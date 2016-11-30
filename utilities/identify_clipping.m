

function identify_clipping(curr_stage,file_range)

    addpath(genpath('./funcs_supporting_local'));
    plot_on = 0;

    if  ~exist('curr_stage','var'), curr_stage = 5; end

    % Path to lfp data
    path_lfp_sample_delay = getpath('path_lfp_sample_delay');
    file_list = get_filelist(path_lfp_sample_delay);
    Nfiles = length(file_list);
    
    if  ~exist('file_range','var')
        file_range = 1:Nfiles;
%         file_range = [60 42];
        %file_range = 43;
    end
    
    
    % Output path
    outpath = fullfile(getpath('path_badchannels'));
    mkdir_dave (outpath);
    outname = ['id_clipping_s' num2str(curr_stage) '.mat'];

    
    parfor i = file_range
        fprintf ('Processing file %d of %d, name %s \n',i,length(file_list),file_list{i});
        [clippingstat{i}, bad_trials{i}, skewout{i}, kurtout{i}, Xint_hist{i}, fnames{i}, lfp_names{i}] = get_clipping_stat(path_lfp_sample_delay,outpath,file_list{i},curr_stage);
    end
    flnames= fname_uname_merge (fnames,lfp_names);  % File and lfp names grouped by file
    
    threshold = get_clipping_threshold_trials;
    [bad_indices, bad_files] = list_badfiles(clippingstat, threshold(2), flnames);
    
    save(fullfile(outpath,outname),'bad_indices','bad_files','flnames','clippingstat', 'bad_trials', 'skewout', 'kurtout', 'fnames', 'lfp_names','threshold');
    
    
    % Merge everything into one long list
    cs_elec= cat(2,clippingstat{:});
    kurt_elec= cat(2,kurtout{:});   
    fln = cat(2,flnames{:});   % File and LFP names for merged list
    %[a b]  = unum2fnum(98,fnames,flnames)          % Restore the original file and lfp numbers from the merged list

    
    if plot_on
        figure; plot(cs_elec); xlabel('Channel'); ylabel('Clipping Statistic');
        figure; hist(cs_elec,500); xlabel('Clipping Statistic'); title('Histogram'); ylabel('N');
    end

end


% Gets a statistic that we will use to classify the file as having
% clipping or not
function [statout, bad_trials, skewout, kurtout, Xint_hist, fname, lfp_names] = get_clipping_stat(inpath,outpath,filename,curr_stage)


    % Note: indices refers to indices of lfp matrix
    % Trials refers to indices of the ctgX_trials arrays
    
    plot_on = 0;
    plot_debug = 0;
    
    if ~exist('inpath','var'); inpath = getpath('path_lfp_sample_delay'); end
    if ~exist('filename','var'); filename = 'L091906.mat'; end
    
    if exist(fullfile(outpath,filename),'file')
        fprintf('File already exists. Skipping.\n');
        return;
    end
    
    % Load Metadata
    md = load (fullfile(getpath('path_metadata'),'sansunits',filename));

    % Load LFP
    load (fullfile(inpath,filename));     % Load full file
    lfp_sample = double(lfp_sample);
    
    % Exctract current stage from metadata
    lfp_stage = get_LFP_ts(lfp_sample,curr_stage,md,lfp_sample_indices);
    Nelect = size(lfp_stage,3);
    dt = get_dt;
    t = (0:size(lfp_stage,1)-1) * dt;
    

    % Instead, I just plot the whole time series and pass it through a
    % filter to take the mean value of the moving window over a fixed time
    % range. Thus, clipping regions summate and average a high value.
    lfpsize = size(lfp_stage);
    intsize = [lfpsize(1)*lfpsize(2) lfpsize(3)];
    Xint = zeros(intsize);
    Xdiff = zeros(intsize);
    
    if nargout>1
        Nhist=100;
        Xint_hist = zeros(Nhist,Nelect);
    end
    
    
    for k = 1:Nelect
        [Xint(:,k), Xdiff(:,k)] = calc_Xint(lfp_stage(:,:,k));
        
        if nargout > 1
            Xint_hist(:,k) = hist(Xint(:,k),Nhist);
        end
    end
    
    lfp_names = md.lfp_names;
    lfp_names=convert_unit_underscores(lfp_names);    % Replace underscores
    skewout = skewness(Xint);
    kurtout = kurtosis(Xdiff);
    fname = filename;
    
    % Prepping some variables
    sz_lfp=size(lfp_stage);
    thresholds = get_clipping_threshold_trials;
    
%     % Test to make sure reshape works like Ithink it does
%     lfp_stage2=reshape(lfp_stage,[sz(1)*sz(2) sz(3)]);
%     lfp_stage3=reshape(lfp_stage2,[sz(1),sz(2),sz(3)]);
%     figure; plot(lfp_stage(:,3,5)); hold on; plot(lfp_stage3(:,3,5),'r.');
    
    bad_trials=[];
    % Find bad trials
    chunk_size_tr = sz_lfp(1);
    [bad_trials] = get_badchunks(Xint,chunk_size_tr);
    Nchunks = size(bad_trials,1);
    num_badtrials = sum(bad_trials);
    percent_badtrials = sum(bad_trials) / Nchunks *100;
    
    % Find bad chunks of sizse 100ms (same as filter length)
    if plot_on  %I was using this as return value, but now only in plots. I've removed it for speed. CAn uncomment if you want to return this later.
        chunk_size = round(0.1/get_dt);
        [bad_chunks] = get_badchunks(Xint,chunk_size);
        Nchunks = size(bad_chunks,1);
        num_badchunks = sum(bad_chunks);
        percent_badchunks = sum(bad_chunks) / Nchunks * 100;
    end
    
    if plot_on
        % Plot raw data and integrated signal
        enum = 6;
        scale_val = std(Xint(:,enum));
        t = 1:size(Xint,1); t=t*dt;
        lfp_stage2 = lfp_stage(:,:,enum); lfp_stage2 = lfp_stage2(:);
        figure; plot(t,zscale(lfp_stage2,scale_val));
        hold on; plot(t,(Xint(:,enum)),'r');
        hold on; plot(t,zscale(2.^-abs(Xdiff(:,enum)),scale_val) - scale_val*10,'k');
        
        
        % Plot chunks data
        Nchunks = size(bad_chunks,1);
        tchunk = get_chunktimes(chunk_size,Nchunks);
        hold on; plot(tchunk,zscale(bad_chunks(:,enum),scale_val),'y','LineWidth',2);
        
        % Plot Trials
        Nchunks = size(bad_trials,1);
        tchunk = get_chunktimes(chunk_size_tr,Nchunks);
        hold on; plot(tchunk,zscale(bad_trials(:,enum),scale_val),'m','LineWidth',2);
        
        
        keyboard

    end
    
    statout = percent_badtrials;
    
    
end



function [bad_chunks] = get_badchunks(Xint,chunk_size)

    Num_electrodes = size(Xint,2);
    Nchunks = floor(length(Xint)/chunk_size);
    Xint2 = Xint(1:Nchunks*chunk_size,:);      % Crop any leftover data on Xint in preparation for chunking
    Xich=reshape((Xint2),[chunk_size, Nchunks, Num_electrodes]);  % Xi chunked
    Xich = squeeze(max(Xich));
    thresholds = get_clipping_threshold_trials;
    bad_chunks = (Xich >= thresholds(1));
    
end


function tchunk = get_chunktimes(chunk_size,Nchunks)
    t=1:chunk_size*Nchunks; t = t*get_dt;
    tchunk = reshape(t,[chunk_size,Nchunks]);
    tchunk = mean(tchunk);
end





