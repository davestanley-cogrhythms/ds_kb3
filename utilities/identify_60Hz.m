


function identify_60Hz(curr_stage,file_range)

    addpath(genpath('./funcs_supporting_local'));
    plot_on = 0;

    if  ~exist('curr_stage','var'), curr_stage = 5; end

    
    path_lfp_sample_delay = getpath('path_lfp_sample_delay');
    file_list = get_filelist(path_lfp_sample_delay);
    Nfiles = length(file_list);
    
    
    if  ~exist('file_range','var')
        file_range = 1:Nfiles;
        %file_range = [42];
    end
    
    
    % Output path
    outpath = fullfile(getpath('path_badchannels'));
    mkdir_dave (outpath);
    outname = ['id_60Hz_s' num2str(curr_stage) '.mat'];

    parfor i = file_range
        fprintf ('Processing file %d of %d, name %s \n',i,length(file_list),file_list{i});
        [linestat{i}, psdout{i}, f{i}, fnames{i}, lfp_names{i}] = get_60Hz_stat(path_lfp_sample_delay,outpath,file_list{i},curr_stage);
    end
    f = f{1};   % Need to do this for parfor
    flnames= fname_uname_merge (fnames,lfp_names);  % File and lfp names grouped by file
    
    threshold = get_60Hz_threshold;
    [bad_indices, bad_files] = list_badfiles(linestat, threshold, flnames);
    
    save(fullfile(outpath,outname),'bad_indices','bad_files','flnames','linestat', 'psdout', 'f', 'fnames', 'lfp_names','threshold')
    
    % Merge everything into one long list
    line_elec= cat(2,linestat{:});
    fln = cat(2,flnames{:});   % File and LFP names for merged list
    
    
    if plot_on
        figure; plot(line_elec); xlabel('Channel'); ylabel('Clipping Statistic');
        figure; hist(line_elec,500); xlabel('Clipping Statistic'); title('Histogram'); ylabel('N');
    end

end






% Gets a statistic that we will use to classify the file as having
% clipping or not
function [lineratio, psdout, f, fname, lfp_names] = get_60Hz_stat(inpath,outpath,filename,curr_stage)


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
    

    % Calculate PSD and take ratio of low freq to 60 Hz. Essentially signal
    % noise ratio
    
    
    % Run test to get dimensionality & initialize output variables
    k=1;
    [temp, psdout_temp, f] = signal_60Hz_ratio_filt (lfp_stage(:,:,k)); clear temp
    psdout = zeros([length(psdout_temp) Nelect]);
    lineratio = zeros(1,Nelect);
    
    % Step through main loop, calculating psd and getting stat ratio
    for k = 1:Nelect
        [lineratio(k), psdout(:,k), f] = signal_60Hz_ratio_filt (lfp_stage(:,:,k));
    end
    
    
    % Define output file names
    lfp_names = md.lfp_names;
    lfp_names=convert_unit_underscores(lfp_names);    % Replace underscores
    fname = filename;
    
end


function [lineratio, P, f] = signal_60Hz_ratio_filt (lfp)

    line_wind = [ 58 62 ] ;
    %sig_wind = [10 50];
    sig_wind = [10 100];

    dt = get_dt;
    N=size(lfp,1);
    %[P f] = pwelch(lfp,N,1,N,1/dt);
    params.tapers = [1 1];  % Use only 1 taper to get the best bandwitdh
    params.trialave = 1;
    [P, f] = psd_wrapper(lfp,'fs',1/dt,'mode',2,'params',params);
    
    
    index60 = f >= line_wind(1) & f <= line_wind(2);
    indexsig = f >= sig_wind(1) & f <= sig_wind(2);
    
    % Take the mean of the two peaks
    lineratio = mean(P(index60)) / mean(P(indexsig));
    
    % Instead, take the fractional power
    lineratio = sum(P(index60)) / sum(P(indexsig));

end
