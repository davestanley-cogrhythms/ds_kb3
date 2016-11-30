

function item_nums = filenames2item_numbers (files_all)
    plot_on = 0;
    
    files_sjefPD = cellfun(@(x) ~isempty(x),strfind(files_all,'sjefPD'));
    files_sgjefPSD = cellfun(@(x) ~isempty(x),strfind(files_all,'sgjefPSD'));
    
    item_nums = zeros(length(files_all),1);
    
    % Process sjefPD files initially
    temp = files_all(files_sjefPD);
    temp = strrep(temp,'sjefPD','');
    temp = strrep(temp,'.bmp','');
    temp = cellfun(@str2num,temp)+1;    % Item number is 1 more than file number. Hence, shift by +1. See item2filename_recording
    item_nums(files_sjefPD) = temp;
    
    % Then process sgjefPSD files
    temp = files_all(files_sgjefPSD);
    temp = strrep(temp,'sgjefPSD','');
    temp = strrep(temp,'.bmp','');
    temp = cellfun(@str2num,temp)+235;    % Item number is 235 more than file number. Hence, shift by +1. See item2filename_recording
    item_nums(files_sgjefPSD) = temp;
    

    if plot_on
        figure; plot(item_nums)
    end
end
