

function adat = func_recalc_pairsmapping_deleteme(sfc_mode,stage_range,file_range)
% Testing scripts - see below

% % Test ulmap - seems to work!
% units_all = horzcat(md.unit_names);
% lfps_all = horzcat(md.lfp_names);
% i = round(unifrnd(1,length(ulpairs)))
% name_Curr = ['Unit:' units_all{ulpairs(1,i)} ' vs LFP:' lfps_all{ulpairs(2,i)}]







%% Function body


 
    % Load SFC for the current stage in stage_range.
    for i = 1:length(stage_range)
        
        % Get file list for specified sfc_mode and stage
        [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
        [fname_suffix] = build_sfcmode(sfc_mode, mode_subgroups);

        mdpath = fullfile(getpath('path_metadata'),'sansunits');
        datapath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(stage_range(i))]);

        if ~exist('file_list','var'); file_list = []; end
        if isempty(file_list); file_list = get_filelist(mdpath); end
        
        if ~exist('file_range','var')
            file_range = [];
        end
        
        if isempty(file_range)
            file_range = 1:length(file_list);
        end

        md = cellfun(@(filename) load (fullfile(mdpath,filename)),file_list);
        %data = cellfun(@load,fullfilec(datapath,file_list(file_range)));
        md=md(file_range);
        
        cum_Ncells = 0;
        cum_Nelects = 0;
        Ncells_shift = [];
        Nelects_shift = [];
        %ulpairs = [];
        %lumap_arr = [];
        
        for j = file_range
            
            % If data file exists, add cumulative unit/electrode count to list of pairs
            if exist(fullfile(datapath,file_list{j}),'file');
                data = load(fullfile(datapath,file_list{j}));
                
                mypairs = data.sfc{1}.lfp_pairs;        % This should be my_pairs for newer data - delete this file and switch to the new version!
                Npairs = size(mypairs,2);

                % Vectors for shifting sfc.mypairs to use absolute (file-independent) numbering scheme
                Ncells_shift = [Ncells_shift cum_Ncells*ones(1,Npairs)];
                Nelects_shift = [Nelects_shift cum_Nelects*ones(1,Npairs)];
            end
            
            
            Ncells = length(md(j).unit_names);
            Nelects = length(md(j).lfp_names);
            
            
            
            % Mapping for units onto electrodes
%             ulpairs = md(j).unit_LFP_pairs;
%             ulpairs2= ulpairs;
%             ulpairs2(1,:) = ulpairs2(1,:) + cum_Ncells;
%             ulpairs2(2,:) = ulpairs2(2,:) + cum_Nelects;
%             ulpairs = [ulpairs ulpairs2];
%             We don't need this anymore, because we have a function that calculates this
            
            % Hash mapping for electrodes on to units (old code)
%             lumap = Inf*ones(1,Nelects);
%             lumap(ulpairs(2,:)) = ulpairs(1,:);
%             lumap = lumap + cum_Ncells;
%             lumap_arr = [lumap_arr lumap];
            
            
            % Update cumulative Ncell and Nelect counts
            cum_Ncells = cum_Ncells + Ncells;
            cum_Nelects = cum_Nelects + Nelects;
            
        end
        
        %Function to calculate and number ulpairs from all files
        ulpairs = build_unit_LFP_pairs_all(md);
        
        % Hash mapping for electrodes on to units (new code)
        Nelects_est = cum_Nelects;
        lumap = Inf*ones(1,Nelects_est);
        lumap(ulpairs(2,:)) = ulpairs(1,:);

        
        map.Ncells_shift = Ncells_shift;
        map.Nelects_shift = Nelects_shift;
        map.ulpairs = ulpairs;
        map.lumap = lumap;

        % Add to struct
        adat.(stagename(stage_range(i))) = map;
    end

    
end


% Old version of function that skips some files
% function adat = func_recalc_pairsmapping(sfc_mode,stage_range,file_range)
%  
%     % Load SFC for the current stage in stage_range.
%     for i = 1:length(stage_range)
%         
%         % Get file list for specified sfc_mode and stage
%         [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
%         [fname_suffix] = build_sfcmode(sfc_mode, mode_subgroups);
% 
%         datapath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(stage_range(i))]);
% 
%         if ~exist('file_list','var'); file_list = []; end
%         if isempty(file_list); file_list = get_filelist(datapath); end
%         
%         if ~exist('file_range','var')
%             file_range = [];
%         end
%         
%         if isempty(file_range)
%             file_range = 1:length(file_list);
%         end
% 
%         md = cellfun(@(filename) load (fullfile(getpath('path_metadata'),'sansunits',filename)),file_list);
%         data_default = cellfun(@load,fullfilec(datapath,file_list(file_range)));
%         
%         cum_Ncells = 0;
%         cum_Nelects = 0;
%         Ncells_shift = [];
%         Nelects_shift = [];
%         ulpairs_arr = [];
%         lumap_arr = [];
%         for j = 1:length(data_default)
%             Ncells = length(md(j).unit_names);
%             Nelects = length(md(j).lfp_names);
%             
%             mypairs = data_default(j).sfc{1}.mypairs;
%             Npairs = size(mypairs,2);
%             
%             % Vectors for shifting sfc.mypairs to use absolute (file-independent) numbering scheme
%             Ncells_shift = [Ncells_shift cum_Ncells*ones(1,Npairs)];
%             Nelects_shift = [Nelects_shift cum_Nelects*ones(1,Npairs)];
%             
%             % Mapping for units onto electrodes
%             ulpairs = md(j).unit_LFP_pairs;
%             ulpairs2= ulpairs;
%             ulpairs2(1,:) = ulpairs2(1,:) + cum_Ncells;
%             ulpairs2(2,:) = ulpairs2(2,:) + cum_Nelects;
%             ulpairs_arr = [ulpairs_arr ulpairs2];
%             
%             % Mapping for electrodes on to units
%             lumap = Inf*ones(1,Nelects);
%             lumap(ulpairs(2,:)) = ulpairs(1,:);
%             lumap = lumap + cum_Ncells;
%             lumap_arr = [lumap_arr lumap];
%             
%             
%             % Update cumulative Ncell and Nelect counts
%             cum_Ncells = cum_Ncells + Ncells;
%             cum_Nelects = cum_Nelects + Nelects;
%             
%         end
%         
%         map.Ncells_shift = Ncells_shift;
%         map.Nelects_shift = Nelects_shift;
%         map.ulpairs_arr = ulpairs_arr;
%         map.lumap_arr = lumap_arr;
% 
%         % Add to struct
%         adat.(stagename(stage_range(i))) = map;
%     end
% 
%     
% end