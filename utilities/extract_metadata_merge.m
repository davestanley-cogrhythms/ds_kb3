
function extract_metadata_merge(metadatapath1,metadatapath2,outpath)
    % Extracts metadata and unit times from Jefferson's new
    % ctgX files and merges them into the original metadata files. Then
    % saves the result in outpath
    
    addpath(genpath('./funcs_supporting_local'));
    
    clc
    data_path = getpath('path_roydata_orig');
    if ~exist('metadatapath1','var'); metadatapath1=fullfile(getpath('path_roydata_orig'),'metadata'); end
    if ~exist('metadatapath2','var'); metadatapath2=fullfile(getpath('path_roy_newctgX')); end
    if ~exist('outpath','var'); outpath=fullfile(getpath('path_metadata')); end

    
    mkdir_dave (outpath);
    mkdir_dave(fullfile(outpath,'sansunits'));
    file_list = dir([fullfile(metadatapath1,'*.mat')]);
    
    % Preload electrode locations data (only works for roy data!)
    if ~get_iscromer
        do_monkeyL = 1; [el_L] = read_in_recording_locations_func(do_monkeyL);
        do_monkeyL = 0; [el_O] = read_in_recording_locations_func(do_monkeyL);
    else
        warning('Electrode location data not being saved. Will be unable to plot electrode locations - not really a big deal.');
    end
    
    for i = 23:23 %length(file_list)
        
        filename = file_list(i).name;
        fprintf ('Processing file %d of %d, name %s \n',i,length(file_list),filename);
        
        s1 = load (fullfile(metadatapath1,filename));
        
        [pathstr, basename, ext] = fileparts(file_list(i).name);
        metadata2_fullname = fullfile(metadatapath2,[basename '_new_conds_trials' ext]);
        clear pathstr versn
        
        
        unmerged = 0;
        if exist(metadata2_fullname,'file')     % Check if new data from Jefferson actually exists!
            
            s2 = load (metadata2_fullname);
            s3 = struct_merge(s1,s2);       % Takes variables in s1, drops them into s2, and saves the output as s3. With overwriting.
            clear s1 s2;
        else
            unmerged=1;
            s3=s1;
            clear s1;
        end
        
        % For cromer data, we need to add 16 to certain electrode numbers
        % This code only executes if the electrode numbers are negative
        s3 = metadata_fix_electrodenumbers(s3);
        
        % Calculate some additional metadata parameters. Converts units
        s3 = metadata_tweaks(s3);  
        
        % Estimate match vs nonmatch trials
        [s3.tr_match,s3.tr_nonmatch] = decode_match_nonmatch(s3);

        % Calculate preferred units. This has to go after metadata_tweaks...
        preferred_path = fullfile(getpath('path_preferred'),'PFC_category_unit_names_CRC.mat');
        if exist(preferred_path,'file')
            [s3.metadata_pref, s3.metadata_prefmat] = calc_metadata_preferred(filename,s3); % Load Jefferson's preferred units into metadata
        else
        end
        
        % Calculate reaction times
        [s3.match_rt, s3.nonmatch_rt] = calc_rt(s3);
        
        % Estimate date of recording from file name
        [s3.datenum_numeric] = est_date(s3);            % Datenum_numeric is in the format output by datenum (i.e. days since 0 AD)
        
        
        % Precalculate morph percentages and schemes
        Ntrials = length(s3.start_times);
        for iii = 1:Ntrials
            [s3.sch{iii}, s3.morphA(iii), s3.morphB(iii)] = condition2morph(s3.each_trial_condition(iii)); % Slow, but not limiting.
        end
        
        % Save electrode locations (only works for Roy data so far!)
        if ~get_iscromer
            if s3.recording_name(1) == 'L'; s3.electrode_locations = extract_electrode_locations(s3,el_L);
            elseif s3.recording_name(1) == 'O'; s3.electrode_locations = extract_electrode_locations(s3,el_O);
            else
                error('Monkey name not found');
            end
        else
        end

        % Save result
        save(fullfile(outpath,filename),'-struct','s3');

        % Save result with units cropped off
        s4 = crop_units_field(s3);         % Adds in neuron preferred schemes (in both struct and mat form)
        save(fullfile(outpath,'sansunits',filename),'-struct','s4');
        
    end
    
end


function s = crop_units_field(s)
    s = rmfield(s,'unit_array'); % Remove unit fields which we don't need
    s = rmfield(s,'unit_ind');   % These take up a lot of space
end




function sout = metadata_tweaks(s)

    fprintf('Running some stuff \n');
    dt = get_dt;
    
    vars_pull(s);
    clear s
    
    % Work in seconds
    [start_times, stop_times, sample_on_times, unit_array, all_encode_times, metadata.units_converted] = convert_to_seconds(start_times, stop_times, sample_on_times, all_encode_times, unit_array);
    
    % Convert seconds to indices
    ss_ind = round([start_times(:) stop_times(:)]/dt);  % % [Start stop] indices
    sample_on_ind = round(sample_on_times / dt);
    for i = 1:length(unit_array)
        unit_ind{i} = round(unit_array{i} / dt);    % % Unit index array
    end
    
    % % Identify electrode / LFP pairs
    [unit_LFP_pairs metadata.unit_LFP_pairs_metadata]= get_unit_LFP_pairs(unit_names,lfp_names);
    
    % % Find trial indices associated with each sample_on_time
    sample_on_trials = get_sample_on_trials(sample_on_ind, ss_ind);
    
    clear i
    sout = vars_push;

    
end


function [start_times, stop_times, sample_on_times, unit_array, all_encode_times, units_converted] = convert_to_seconds(start_times, stop_times, sample_on_times, all_encode_times, unit_array)

    % This function checks if the units of start_times, stop_times,
    % sample_on_times, and unit_array are in milliseconds. If they are,
    % convert it to seconds

    % The difference between the start and stop times should be about 
    % 10 seconds (certainly less than 500)
    units_converted.readme = 'Unit conversions that were made';
    if median(stop_times-start_times) > 500
        % Work in milliseconds
        start_times = start_times / 1e3;
        stop_times = stop_times / 1e3;
        units_converted.start_stop_times = 1;
    else
        units_converted.start_stop_times = 0;
    end
    
    % The difference between subsequent sample_on_times should be about 
    % 10 seconds (certainly less than 500)
    if median(diff(sample_on_times)) > 500
        % Work in milliseconds
        sample_on_times = sample_on_times / 1e3;
        units_converted.sample_on = 1;
    else
        units_converted.sample_on = 0;
    end
    
    % The mean interspike interval should be about 0.1 seconds (certainly
    % less than 5 seconds)
    if mean(diff(unit_array{1})) > 5
        % Work in milliseconds
        for i = 1:length(unit_array)
            unit_array{i} = unit_array{i} / 1e3;
        end
        units_converted.unit_array = 1;
    else
        units_converted.unit_array = 0;
    end
    
    if median(diff(all_encode_times)) > 0.1         % Spacing between many event codes is only a few milliseconds!
        all_encode_times = all_encode_times / 1e3;
        units_converted.all_encode_times = 1;
    else
        units_converted.all_encode_times = 0;
    end
end


function s3 = metadata_fix_electrodenumbers(s3)
    lfp_names = s3.lfp_names;
    temp = strrep(lfp_names,'PFC_Elect_','');
    lfp_nums = cellfun(@str2num,temp);
    
    if any(lfp_nums < 1)
        lfp_nums = lfp_nums + 16;

        for i = 1:length(lfp_names)
            lfp_names2{i} = ['PFC_Elect_' num2str(lfp_nums(i))];
        end

        s3.lfp_names = lfp_names2;
    end
    
end




function [match_rt, nonmatch_rt] = calc_rt(md)


    % % Define constants

    roys_way = 0;

    if ~get_iscromer
        code_shift = 3;
        code_match_response_shift = 3;
        code_nonmatch_response_shift = 5; % I think
        code_sample_on = 27;
        code_test = 29;
        code_short = 4;
        code_long = 30;
    else
    %     % Old method ... this one catches some bugs
    %     code_shift = 3;
    %     code_sample_on = 23;
    %     code_test = 25;
    %     code_short = 175;
    %     code_long = 26;

        code_shift = 4;
        code_match_response_shift = 3;
        code_nonmatch_response_shift = 5;
        code_sample_on = 23;
        code_test = 25;
        code_short = 176;
        code_long = 27;
    end



    % % Run analysis

    clear match_rt nonmatch_rt schA schB schA2
    
    vars_pull(md);


    % Find all instances of 27 - these are sample on indices
    sample_on_encodes=find(all_encodes == code_sample_on);


    %Error trials
    etrs = each_trial_response(sample_on_trials);
    ind_err = etrs ~= 0;     % Trials for which the correct answer was not given (for whatever reason)

    if roys_way
        inds = logical(mod(each_trial_condition,2));        % Odd for match
        indl = ~logical(mod(each_trial_condition,2));       % Even for non-match
        inds = inds(sample_on_trials) & ~ind_err;
        indl = indl(sample_on_trials) & ~ind_err;
    else % My way, using all_encodes
        encodesp4 = all_encodes(sample_on_encodes+code_shift);
        inds = encodesp4 == code_short & ~ind_err;
        indl = encodesp4 == code_long & ~ind_err;
    end

    ind_rest = ~(inds | indl | ind_err);

    %match_dt = all_encode_times(ind2+code_match_response_shift-1)-all_encode_times(ind2); % test time since sample on - should be about 1.6 seconds
    match_rt = all_encode_times(sample_on_encodes+code_match_response_shift)-all_encode_times(sample_on_encodes+code_match_response_shift-1);
    match_rt(~inds | ind_err) = NaN;

    %nonmatch_dt = all_encode_times(ind2+code_nonmatch_response_shift-1)-all_encode_times(ind2); % Test time since sample on - should be about 3.2s
    nonmatch_rt = all_encode_times(sample_on_encodes+code_nonmatch_response_shift)-all_encode_times(sample_on_encodes+code_nonmatch_response_shift-1);
    nonmatch_rt(~indl | ind_err) = NaN;



end


function [date_out] = est_date(md)

    recording_name = md.recording_name;
    if ~get_iscromer
        date_string = recording_name(2:7);  % Jeffersons format is mm/dd/yy
        
        mm = str2double(date_string(1:2));
        dd = str2double(date_string(3:4));
        yy = str2double(date_string(5:6));
        
    else
        date_string = recording_name(3:8);  % Cromer's format is yy/mm/dd
        yy = str2double(date_string(1:2));
        mm = str2double(date_string(3:4));
        dd = str2double(date_string(5:6));
        
    end
    
    date_out = datenum(yy+2000,mm,dd);
    

end



function electrode_locations = extract_electrode_locations(s3,el_struct)

    if ~get_iscromer
        date_string = s3.recording_name(2:7);  % Jeffersons format is mm/dd/yy
    else
        date_string = s3.recording_name(3:8);  % Cromer's format is yy/mm/dd
    end
    
    % Find matching dates
    mydate = str2num(date_string);
    inds = find(el_struct.dateinfo == mydate);
    
    for i = 1:length(s3.lfp_names)  % For each electrode...
        temp = s3.lfp_names{i};
        lfp_num = str2num(strrep(temp,'PFC_Elect_',''));
        inds2 = el_struct.electrode_num(inds) == lfp_num;               % Find corresponding electrode in el_struct
        electrode_locations.AP(i) = el_struct.AP(inds(inds2));          % Save corresponding coordinates ...
        electrode_locations.ML(i) = el_struct.ML(inds(inds2));
        electrode_locations.electrode_num(i) = el_struct.electrode_num(inds(inds2));    % ... and electrode number.
        
    end
    
end


