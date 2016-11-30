

function ind = filename2cond_index(cond_filenames,filename)
    cond_filenames = strrep(cond_filenames,'cond_o','');
    cond_filenames = strrep(cond_filenames,'cond_l','');
    cond_filenames = strrep(cond_filenames,'.con','');
    filename = strrep(filename,'L','');
    filename = strrep(filename,'O','');
    filename = strrep(filename,'.mat','');
    
    cond_num = cellfunu(@str2num,cond_filenames);
    cond_num = cell2mat(cond_num);
    
    filename_num = str2num(filename);
    
    ind = find(filename_num == cond_num);
    if isempty(ind)
        warning(['Corresponding conditions file not found for file ' filename ' !']);
        ind = 2;    % Temporarily use index = 2 for missing files. Jefferson hopefully can send me this
    end
end