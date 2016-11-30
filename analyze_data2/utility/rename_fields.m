


% Fix missing files
% Add missing fields to structures as needed

%% Setup vars
sfc_mode=3.20001;
curr_stage=5;
field_to_rename = {'S2ave','S2AVE'};

%% First, get list of files to fix

[~, mode_subgroups] = decode_sfc_mode(sfc_mode);
[fname_suffix] = build_sfcmode(sfc_mode, mode_subgroups);
outpath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(curr_stage)]);


fn = dir(fullfile(outpath,'*.mat'));
fn = {fn.name};


%%
for f = 1:length(fn)
    
    
    temp = load(fullfile(outpath,fn{f}));
    sfc = temp.sfc;                     % Necessary for running with parfor - needs to know sfc is a variable.
    changes_made=0;
    for i = 1:length(sfc)
        if ~isfield(sfc{i},field_to_rename{2})
            fprintf(['Field ' field_to_rename{2} ' is missing from file # ' num2str(f) ',' fn{f} ', condition j=' num2str(i) '\n']);
            sfc{i}.(field_to_rename{2}) = sfc{i}.(field_to_rename{1});
            sfc{i} = rmfield(sfc{i},field_to_rename{1});
            
%             %Tweak some additional fields as well ... should comment it out if just want to rename 1 field
%             sfc{i}.('ConfC') = sfc{i}.('confC');
%             sfc{i} = rmfield(sfc{i},'confC');
%             sz=size(sfc{i}.Ssp);
%             sfc{i}.ConfC = reshape(sfc{i}.ConfC,[1, sz(2:3)]);
%             changes_made=1;
        end
    end
    
    if changes_made
        save(fullfile(outpath,fn{f}),'sfc');        % Parfor doesn't work with save :(
    end
    clear sfc
    
end
%%

    



