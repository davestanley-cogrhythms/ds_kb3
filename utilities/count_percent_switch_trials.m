

function count_percent_switch_trials (file_range)

% % % Calculates the percent switch trials and saves to a file.
% I have another function to do this now, called in plots_preferred, so probably don't need this
% anymore.
% % % 
format compact
format long g

run_setdefaultfig

plot_debug = 1;


% path_metadata = getpath('path_metadata');
path_lfp_sample_delay = getpath('path_lfp_sample_delay');
% path_lfp = getpath('path_lfp');


file_list = get_filelist(path_lfp_sample_delay);
Nfiles = length(file_list);


if  ~exist('file_range','var')
    file_range = 1:Nfiles;
%     file_range = [25:30];
%     file_range = [40];
end

for i = file_range
    fprintf ('Processing file %d of %d, name %s \n',i,length(file_list),file_list{i});
    percent_switches(i) = calc_stats(file_list{i});
end

outfile = fullfile(getpath('path_buffer_curr'),'percent_switches');
save(outfile,'percent_switches','file_range','file_list');

if plot_debug
    figure; plot(file_range,percent_switches);
    xlabel('File number'); ylabel('Percent Switch Trials')
end

end




function percent_switches = calc_stats(filename)

    if ~exist('filename','var'); filename = 'L091906.mat'; end
    
    
    include_switch_trials = 1;      % Include in ctgsetli 
    
    md = load (fullfile(getpath('path_metadata'),'sansunits',filename));     % Load metadata
    %load (fullfile(getpath('path_lfp'),filename));     % Load full file
    %lfpmat = double(lfpmat_int16); clear lfpmat_int16;
    
    
    Ntrials = length(md.start_times);  % Total number of trials
  
    if include_switch_trials == 1
        clear sch
        for iii = 1:Ntrials
            [sch(iii)] = condition2morph(md.each_trial_condition(iii)); % Slow, but not limiting.
        end
        
        sch = (sch == 'A');
        sch_curr = sch(2:end);
        sch_prev = sch(1:end-1);

        tsw = sch_curr ~= sch_prev; % Trials switch
        tct = sch_curr == sch_prev; % Trials continuous

        tsw = logical([0 tsw]);
        tct = logical([0 tct]); % First trial is neither switch nor continuous - they're nothing!
        
        % Reduce switch trials to be indexed by samples
        tsw = tsw(md.sample_on_trials);
        tct = tct(md.sample_on_trials);
        
        % Display percent switch trials
        Ntsw = sum(tsw); Ntct = sum(tct);
        fprintf('Percent switch trials = %g \n', (Ntsw / (Ntsw + Ntct) *100) );
        percent_switches = (Ntsw / (Ntsw + Ntct) *100);
    end
    
    

end




