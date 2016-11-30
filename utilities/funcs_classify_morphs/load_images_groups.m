

function [A, file_list] = load_images_groups(path)
%% out = load_images_groups(path)

    file_list = get_filelist(path,'*.bmp');
    
    for i = 1:length(file_list)
        A{i} = imread(fullfile(path,file_list{i}));
%         pause
    end
    
    

end