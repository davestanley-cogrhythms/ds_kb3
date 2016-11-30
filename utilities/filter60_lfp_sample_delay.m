
function filter60_lfp_sample_delay(file_range)

    inpath = fullfile(getpath('path_lfp_sample_delay'));
    outpath = fullfile(getpath('path_lfp_sample_delay'),'filtered_temp');
    
    mkdir_dave(outpath);
    
    file_list = dir(fullfile(inpath,'*.mat'));
    Nfiles = length(file_list);

    if  ~exist('file_range','var')
        file_range = [1:Nfiles];
    end
    
    for i = file_range
        fprintf ('Processing file %d of %d, name %s \n',i,length(file_list),file_list(i).name);
        filter_file(file_list(i).name,inpath,outpath);
    end
    
end


function filter_file(filename, inpath, outpath)


    % Note: indices refers to indices of lfp matrix
    % Trials refers to indices of the ctgX_trials arrays
    
    
    if ~exist('filename','var'); filename = 'L091906.mat'; end
    
    
    
    if exist(fullfile(outpath,filename),'file')
        fprintf('File already exists. Skipping.\n');
        return;
    end

    load (fullfile(inpath,filename));     % Load unfiltered data
    lfp_sample = double(lfp_sample);
    lfp_sample_orig = lfp_sample;
    
    [lfp_sample, power_reduction, interval_lengths] = smartfilter(lfp_sample, 1/get_dt);
    
    % Format of lfp_sample is lfp_sample(sample data, sample number, electrode number)
    lfp_sample = int16(lfp_sample);
    save(fullfile(outpath,filename),'lfp_sample','metadata_lfp','lfp_sample_indices');
    
    mkdir_dave(fullfile(outpath,'filter_metadata'));
    save(fullfile(outpath,'filter_metadata',filename),'power_reduction','interval_lengths')
    
    
 
end



function [xfilt, power_reduction, interval_lengths] = smartfilter(x, fs)

    plot_on = 0;
    dt = 1/fs;
    
    t = 1:size(x,1); t=t/fs;
    sz = size(x);
    
    tic
    interval_lengths = zeros(sz(2),sz(3));
    xfilt = zeros(sz);
    for k = 1:sz(3)
        for j=1:sz(2);
            [intervals] = smartfilter_interval (t, x(:,j,k));
            [xfilt(:,j,k)] = qif(t,x(:,j,k),intervals);
            
            interval_lengths(j,k) = size(intervals,1);
        end
    end
    toc
    xfilt = floor(xfilt);
    
    varpre = squeeze(var(x));
    varpost = squeeze(var(xfilt));
    power_reduction = varpost ./varpre * 100;
    
    if plot_on
        %% Plotting
        Ntraces = 10;
        traces = floor(unifrnd(1,sz(2),1,Ntraces));
        traces = sort(traces);
        enum=7;
        
        
        figure; 
        for j = 1:length(traces)
            clf;
            fprintf(['Plotting ' num2str(j) 'st trace, ' num2str(traces(j)) '.\n']);
            xcurr = x(:,traces(j),enum);
            xfiltcurr = xfilt(:,traces(j),enum);
            
            N = length(xcurr);
            [pxx, f] = pwelch(xcurr,[round(N)],[1],[round(N)],1/dt);
            subplot(211); plot(t,xcurr,'b');
            subplot(212); semilogy(f,pxx,'b'); 

            [pxx_filt, f_filt] = pwelch(xfiltcurr,[round(N)],[1],[round(N)],1/dt);
            subplot(211); hold on; plot(t,xfiltcurr,'r');
            subplot(212); hold on; plot(f,pxx_filt,'r'); 
            
            legend('Original',['Filt ' num2str(power_reduction(traces(j),enum)) '% powr orig']); 
            pause;
        end
    end

end