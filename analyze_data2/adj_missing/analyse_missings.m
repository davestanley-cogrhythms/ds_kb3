
stage = 3;
mode = 2.2;

adj_mode = '_adj';
datapath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(mode) adj_mode],['stage_' num2str(stage)]);

file_list = get_filelist(datapath); 
file_range = 1:79;

% Load sfc data
data_default = cellfun(@load,fullfilec(datapath,file_list(file_range)));

% Load md data
file_list = get_filelist(fullfile(getpath('path_metadata'),'sansunits'));
md = cellfun(@(filename) load (fullfile(getpath('path_metadata'),'sansunits',filename)),file_list(file_range));
    

% Extract missing data
for i = 1:79
    file_missing{i} = data_default(i).sfc{1}.adj_was_missing;
    nmissings(i) = sum(file_missing{i});
end

% Count number of missing units per file
fnums_missing = find(nmissings > 0);
fnames_missing = file_list(fnums_missing);

for i = 1:length(fnums_missing)
    curr_fnum = fnums_missing(i);
    units_missing{i} = md(curr_fnum).unit_names(logical(file_missing{curr_fnum}));
end

% % Extract missing unit names
% fu_missing = logical(cat(2,file_missing{:}));
% funames = cat(2,md(:).unit_names);
% 
% for i = 1:length(fnums_missing)
%     unit_names_missing{i} = funames(fu_missing)
% end



