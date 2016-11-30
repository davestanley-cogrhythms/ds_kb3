

% sfc_mode = 22.40171110; merge_i0(sfc_mode,3);
% sfc_mode = 52.70030100; merge_i0(sfc_mode,4);


function merge_i0(sfc_mode,curr_stage)


    %% Setup vars
    if ~exist('sfc_mode','var'); sfc_mode=22.4415111; end
    if ~exist('curr_stage','var'); curr_stage=3; end

    %% First, get list of files to fix

    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix] = build_sfcmode(sfc_mode, mode_subgroups);
    outpath_part1 = fullfile(getpath('path_buffer_curr'),['test4mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(curr_stage)]);
    outpath_part2 = fullfile(getpath('path_buffer_curr'),['test3mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(curr_stage)]);
    outpath_saving = fullfile(getpath('path_buffer_curr'),['test5mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(curr_stage)]);
    mkdir_dave (outpath_saving);


    fn = dir(fullfile(outpath_part1,'*.mat'));
    fn = {fn.name};


    %% Perform the operation
    parfor f = 1:length(fn)
        myconvert(fn,f,outpath_part1,outpath_part2,outpath_saving)
    end
    
end


function myconvert(fn,f,outpath_part1,outpath_part2,outpath_saving)

    

    clear sfc
    fprintf(['Loading file ' fullfile(outpath_part1,fn{f}) '\n'])
    part1 = load(fullfile(outpath_part1,fn{f}));
    part2 = load(fullfile(outpath_part2,fn{f}));
    sfc1=part1.sfc;
    sfc2=part2.sfc;
    sfc = sfc2;

    for i = 1:length(sfc)

        if i <= length(sfc1)
            sfc{i} = sfc1{i};
        end
    end
    
    save(fullfile(outpath_saving,fn{f}),'sfc');
end