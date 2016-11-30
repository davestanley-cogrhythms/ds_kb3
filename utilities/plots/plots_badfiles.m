
    % Plots bad clipping and bad 60 Hz data
    % Possibly can delete this file, as it's outdated by
    % plots_badfiles_bytrials.m

    %% Prep data and load (slow!!)
    plot_debug =1;
    plot_on = 1;

    addpath(genpath('./funcs_supporting_local'));
    addpath(genpath('../funcs_supporting_local'));
    
    % Path to lfp data
    path_lfp_sample_delay = getpath('path_lfp_sample_delay');
    file_list = get_filelist(path_lfp_sample_delay);
    Nfiles = length(file_list);

    %Script general parameters
    plot_on = 0;
    curr_stage = 5;
    file_range = 1:Nfiles;
    
    % Path to bad channels for both clipping and 60 Hz
    path_badchannels = fullfile(getpath('path_badchannels'));
    name_clip = ['id_clipping_s' num2str(curr_stage) '.mat'];
    name_60Hz = ['id_60Hz_s' num2str(curr_stage)  '.mat'];
    
    
    % Load bad data info
    s_clip = load(fullfile(path_badchannels,name_clip));
    s_60Hz = load(fullfile(path_badchannels,name_60Hz));
 
    % Merge files
    li_clip = cat(2,s_clip.clippingstat{:});
    li_60 = cat(2,s_60Hz.linestat{:});
%     li_hist = cat(2,s_clip.Xint_hist{:});         % No longer saved
    li_psd = cat(2,s_60Hz.psdout{:});
    li_f = s_60Hz.f;
    
    % % Import filenames (optional)
    fnames = s_clip.fnames;
    flnames = s_clip.flnames;
    [a b]  = unum2fnum(100,fnames,flnames);          % Restore the original file and lfp numbers from the merged list
    

    
    % % Load all LFP data into memory
    % Load Metadata
    md = cellfun(@(filename) load (fullfile(getpath('path_metadata'),'sansunits',filename)),file_list(file_range));

    % Load LFP and merge files (This is HUGE!)
    lfp = cellfun(@(filename) load (fullfile(path_lfp_sample_delay,filename)), file_list(file_range)); % Load all files. This could take a while... this is huge!
    lfp_all = merge_allfiles_lfp(lfp); clear lfp
    
    % Reshape lfp into single rows, and then conver to a matrix. Note, we are cropping the matrix to the length of the shortest channel to ensure that dims are consistent
    lfp_mat = lfp_cells_to_matrix(lfp_all);
    lfp_mat = single(lfp_mat);
    dt = get_dt;
    t = 0:size(lfp_mat,1)-1; t=t*dt;
    Nelect = size(lfp_mat,2);
    
    % % % % Finished loading stuff. Now Make some plots. % % % % %
    
    
    % Prep clipping and 60 Hz data
    thresh_clip = get_clipping_threshold;
    thresh_60Hz = get_60Hz_threshold;
    
    %thresh_clip = 2;  % (For testing only!!!)
    %thresh_60Hz = 0.4;
    
    bad_clip = li_clip > thresh_clip;
    bad_60 = li_60 > thresh_60Hz;
    
    percent_bad_clipping = sum(bad_clip)/length(bad_clip)*100
    percent_bad_60 = sum(bad_60)/length(bad_60)*100
    percent_both = sum(bad_60 | bad_clip) / length(bad_clip) * 100
    
    %% Plot clipping and 60 Hz statistics
    if plot_debug
        figure; 
        subplotsq(4,1); plot(li_clip); title(['Percent bad trials due to clipping ' num2str(percent_bad_clipping)]);
        subplotsq(4,2); hist(li_clip,500)
        
        li_60_temp =  li_60;
        li_60_temp (li_60_temp > 100) = 100; %Set a saturation point so the plot doesn't look stupid
        subplotsq(4,3); plot(li_60_temp); title(['Bad percentage due to 60 Hz ' num2str(percent_bad_60)]);
        subplotsq(4,4); hist(li_60_temp,500)
    end

    %% Plot traces with 60Hz contamination in a certain range
    testrange = [ 0.38 0.42] ; % I am using 40% as my upper limit for 60 Hz contamination. (I.e. at most 40% of 10-100Hz power can result from line noise).
    testind = li_60 >= testrange(1) & li_60 <= testrange(2);
    plot_range = 1:2500;
    if plot_on
        figure; comments = plott_ani_pairs(lfp_mat(plot_range,testind),@(X) filtplot(t(plot_range),X),li_psd(:,testind),@(X) plot(li_f,X),'randcol',1);
    end
    
    %% Plot traces with clipping contamination in a certain range
    testrange = [ 1.3 1.5] ; 
    testind = li_clip >= testrange(1) & li_clip <= testrange(2) & ~bad_60;     % Don't count files that are already flagged as bad for 60 Hz
    %testind = li_clip >= testrange(1) & li_clip <= testrange(2);
    plot_range = 1:1320300;
    if plot_on

        figure;
        %comments = plott_ani_pairs(lfp_mat(plot_range,testind),@(X) intplot(t,X),li_psd(:,testind),@(X) plot(li_f,X));  %Plot with power spectrum
        %comments = plott_ani_pairs(lfp_mat(plot_range,testind),@(X) plot(t(plot_range),X));   % Just the trace
        comments = plott_ani_pairs(lfp_mat(plot_range,testind),@(X) intplot(t(plot_range),X)); xlabel('t(s)');ylabel('LFP');   % As above, but without PSD
        
        
        % Make my own plott_ani to provide more details
        
        testlist = find(testind);
        testlist = testlist(testlist >= 293);
        figure;
        for i = 1:length(testlist);
            clf
            [a b]  = unum2fnum(testlist(i),fnames,flnames);
            subplotrows(2,1); intplot(t(plot_range),lfp_mat(plot_range,testlist(i))); title(['File #' num2str(a) 'e' num2str(b) ' ' fnames{a} ' Skewness = ' num2str(li_clip(testlist(i)))]);
%             subplotrows(2,2); bar(li_hist(:,testlist(i)));
            
            prompt = 'Type a comment or hit enter to continue q to quit: ';
            returns{i} = input(prompt,'s');
            if strcmp(returns{i},'q') || strcmp(returns{i},'Q'); break; end
        end
    end
    
    
    
    
    