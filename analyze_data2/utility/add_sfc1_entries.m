


% When I want to run calculations only on sfc{11}, code further down the
% pipeline will cholk. Therefore, I need to add dummy fields to sfc{1}.
% This is a cheap hack that should go away once I restructure run_analysis
% to produce a "firsts" field that will store all run's metadata.

% 
% This essentially copies all fields from source_cell over to dest_cell in
% all files.

%% Setup vars
sfc_mode=22.401310;
curr_stage=2;
source_cell = 11;
dest_cell = 1;

%% First, get list of files to fix

[~, mode_subgroups] = decode_sfc_mode(sfc_mode);
[fname_suffix] = build_sfcmode(sfc_mode, mode_subgroups);
outpath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(curr_stage)]);


fn = dir(fullfile(outpath,'*.mat'));
fn = {fn.name};


%%
for f = 1:length(fn)
    
    clear sfc
    load(fullfile(outpath,fn{f}));
    
    sfc_src = sfc{source_cell};
    fieldn = fieldnames(sfc_src);
    for i = 1:length(fieldn)
        sfc{dest_cell}.(fieldn{i}) = sfc_src.(fieldn{i})*NaN;
    end
    
    save(fullfile(outpath,fn{f}),'sfc');
end

%%
    
    



