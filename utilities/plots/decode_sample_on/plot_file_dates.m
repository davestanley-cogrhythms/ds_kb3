
%% Dates of all files, extracted from metadata
why
why
why

clear
% close all

addpath(genpath('../../funcs_supporting_local/'));

path_md = fullfile(getpath('path_metadata'),'sansunits');
%path_md = fullfile(getpath('path_metadata'));
%path_md = fullfile(getpath('path_roydata_orig'),'metadata');
%path_md = getpath('path_roy_newctgX');

file_list = get_filelist(path_md);

md_all = cellfun(@load,fullfile(path_md,file_list));


%% Extract dates

dates_all = []; datesnum_all=[];
for i = 1:length(md_all)
    md = md_all(i);
    
    
    recording_name = md.recording_name;
    if ~get_iscromer
        date_string = recording_name(2:7);  % Jeffersons format is mm/dd/yy
        date_string_from_filename = file_list{i}(2:7);
        if ~strcmp(date_string,date_string_from_filename)
            warning('Date in md.recording_name does not match date in file name');
            keyboard
        end
        
        mm = str2double(date_string(1:2));
        dd = str2double(date_string(3:4));
        yy = str2double(date_string(5:6));
        
        
        
    else
        date_string = recording_name(3:8);  % Cromer's format is yy/mm/dd
        date_string_from_filename = file_list{i}(3:8);
        if ~strcmp(date_string,date_string_from_filename)
            warning('Date in md.recording_name does not match date in file name');
            keyboard;
        end
        
        yy = str2double(date_string(1:2));
        mm = str2double(date_string(3:4));
        dd = str2double(date_string(5:6));
        
    end
    
    dates_all = [dates_all [yy; mm; dd]];
    
    
    datesnum_all = [datesnum_all datenum(yy+2000,mm,dd)];
    
end

%% Plot Results

figure; plot(datesnum_all - min(datesnum_all))

