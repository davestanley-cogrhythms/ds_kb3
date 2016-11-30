
% Renames 
file_range = 1:79
for curr_stage = -2:5
    path_badchannels = fullfile(getpath('path_badchannels'));
    name_clip2 = ['id_clipping_s' num2str(curr_stage) '.mat'];
    name_60Hz2 = ['id_60Hz_s' num2str(curr_stage)  '.mat'];


    name_clip = ['id_clipping_s' num2str(curr_stage) '_files' num2str(min(file_range)) '-' num2str(max(file_range)) '.mat'];
    name_60Hz = ['id_60Hz_s' num2str(curr_stage) '_files' num2str(min(file_range)) '-' num2str(max(file_range)) '.mat'];
    
    
    
    %fprintf(['!mv' fullfile(path_badchannels,name_clip) ' ' fullfile(path_badchannels,name_clip2) '\n']);
    eval(['!mv ' fullfile(path_badchannels,name_clip) ' ' fullfile(path_badchannels,name_clip2)]);

end
    