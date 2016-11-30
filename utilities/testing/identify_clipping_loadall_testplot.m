
function identify_clipping_loadall_testplot(curr_stage,file_range)

    addpath(genpath('../funcs_supporting_local'));

    if  ~exist('curr_stage','var'), curr_stage = 3; end

    outpath = fullfile(getpath('path_clipping'));
    mkdir_dave (outpath);
    
    path_lfp_sample_delay = getpath('path_lfp_sample_delay');
    file_list = get_filelist(path_lfp_sample_delay);
    Nfiles = length(file_list);

    if  ~exist('file_range','var')
        file_range = 1:Nfiles;
        file_range = [60];
    end
    
    
    for i = file_range
        fprintf ('Processing file %d of %d, name %s \n',i,length(file_list),file_list{i});
        lfp_diffhist = get_diffhist(path_lfp_sample_delay,outpath,file_list{i},curr_stage);
    end

end


function lfp_diffhist = get_diffhist(inpath,outpath,filename,curr_stage)


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
    md = load (fullfile(getpath('path_metadata'),filename));

    % Load LFP
    load (fullfile(getpath('path_metadata'),filename));     % Load metadata
    load (fullfile(inpath,filename));     % Load full file
    lfp_sample = double(lfp_sample);
    
    % Exctract current stage from metadata
    lfp_stage = get_LFP_ts(lfp_sample,curr_stage,md,lfp_sample_indices);
    t = (0:size(lfp_stage,1)-1) * dt;
    
    
    % Basic test - just take STD of dataset
    Xin = lfp_stage(:,1:4,1:2);
    [tout Xout] = calc_ts_func(t,Xin,0.05,0.9,0.9,@std,0);
    
    
    % Try applying historgram to 200ms segments of the whole dataset.
    % Didn't work out too well, because dataset size is too small.
    [N Ntrials, Nelect] = size(lfp_stage);
    [tout Xoutbins] = calc_ts_func(t,lfp_stage(:,:,:),0.2,0.5,0.9,@(X) histc(diff(X),linspace(min(X(:)),max(X(:)),20)),1);
    [tout Xout] = calc_ts_func(t,lfp_stage(:,:,:),0.2,0.5,0.9,@(X) X,1);
    
    Xout = cellfunu(@(X) permute(X,[4 1 2 3]),Xout);
    Xoutbins = cellfunu(@(X) permute(X,[4 1 2 3]),Xoutbins);
    
    Xout=cell2mat(Xout);
    Xoutbins=cell2mat(Xoutbins);
    Xout=permute(Xout,[2,1,3,4]);
    Xoutbins=permute(Xoutbins,[2,1,3,4]);
    ind1=1:1;ind2=1;
    if plot_debug
        for j=ind1
            figure(1);plott_ani_pairs(Xout(:,:,j,ind2),@plot,Xout(:,:,j,ind2),@(X) hist(diff(X)),Xoutbins(:,:,j,ind2),@plot);    % Cycle through the subsets of each trial
            figure(2);plot(lfp_stage(:,j,ind2)); title('Entire trial');
        end
    end
     
    
    
    % Just calculates histogram independently for each trial, and then
    % plots them. This looks much better - non gaussian plots could be our
    % traces with clipping
    Xe1=lfp_stage(:,:,1);
    Xe1diff = diff(Xe1,1);
    feature=hist(Xe1diff,100);
    feature=arrayfunu( @(dim) hist(Xe1diff(:,dim),100)',1:size(Xe1diff,2));
    feature=cell2mat(feature);
    if plot_on
        figure(1);plott_ani_pairs(feature,@(X) plot(X,'b.','Marker','+','LineStyle','-'),Xe1,@plot, Xe1diff,@plot)
    end

    
    % Here I reshape the matrix into smaller ~200ms segments and plot my nonlinear filter
    % It doesn't really work too well, because there is a lot of noise.
    % Every point of inflection triggers a false alarm.
    Xe1=lfp_stage(:,:,1);
    seglength = 200;
    Xe1=Xe1(:); Xe1=reshape(Xe1,seglength,length(Xe1)/seglength);
    Xfilt = filt_amplify_linearities(Xe1);      % Diff followed by transform X -> 2^-abs(X)
    if plot_debug;
        figure(1);plott_ani_pairs(Xe1,@plot,diff(Xe1,2),@plot,Xfilt,@plot);
    end
    
    % Instead, I just plot the whole time series and pass it through a
    % filter to take the mean value of the moving window over a fixed time
    % range. Thus, clipping regions summate and average a high value.
    time_wind = 0.100; % in ms
    Xe1=lfp_stage(:,:,1);
    Xfilt = filt_amplify_linearities(Xe1);      % Diff followed by transform X -> 2^-abs(X)
    Xint = integrate_sig(Xfilt(:),round(time_wind/dt),dt);
    Xdiff = diff_pad(Xe1,2,1,0);
    
    if plot_on
        figure; h0 = plot(zscore(Xe1(:)));
        hold on ; h1 = plot(zscore(Xint(:)),'r','LineWidth',2)
        hold on ; h2 = plot(zscore(Xfilt(:))-10,'k','LineWidth',1)
        hold on ; h3 = plot(zscore(Xdiff(:))+10,'g','LineWidth',1)
        figure; hist(Xint,100)
    end
    skewness(Xint)
    
end


