

% Merges _i0j0 files into single mat file

%% Example code
% merge_i0j0(22.40131,2); merge_i0j0(22.40131,3); 
% merge_i0j0(2.20111,2); merge_i0j0(2.20111,3); 
% merge_i0j0(2.3015111,2);
% merge_i0j0(2.3015111,3); 
% merge_i0j0(2.20151,2); merge_i0j0(2.20151,3); 
% merge_i0j0(2.201511,2); %merge_i0j0(2.201511,3); 
% merge_i0j0(2.201112,2); %merge_i0j0(2.201112,3); 
% merge_i0j0(22.401311,2); %merge_i0j0(22.401311,3); 
% merge_i0j0(41.621310,2); merge_i0j0(41.621310,3); 
% merge_i0j0(2.221510,2); merge_i0j0(2.221510,3); 
% merge_i0j0(22.421310,2); merge_i0j0(22.421310,3); 
% merge_i0j0(22.401310,2); merge_i0j0(22.401310,3); 
% merge_i0j0(2.201112,2); merge_i0j0(2.201112,3); 
% merge_i0j0(2.201511,2); merge_i0j0(2.201511,3);
% merge_i0j0(42.60121,3);
% merge_i0j0(42.601311,3);
% merge_i0j0(43.60121,3);
% merge_i0j0(43.601311,3);
% merge_i0j0(22.441310,3); merge_i0j0(22.441310,2); 
% merge_i0j0(22.411310,3); merge_i0j0(22.411310,2); 
% merge_i0j0(22.441310,3); merge_i0j0(22.441310,2); 
% merge_i0j0(41.641310,3); merge_i0j0(41.641310,2); 
% merge_i0j0(41.611310,3); merge_i0j0(41.611310,2); 
% merge_i0j0(41.641310,3); merge_i0j0(41.641310,2); 
% merge_i0j0(22.441311,3); merge_i0j0(22.441311,2); 
% merge_i0j0(41.641311,3); merge_i0j0(41.641311,2); 
% merge_i0j0(22.4413101,3); merge_i0j0(22.4413101,2); 
% merge_i0j0(22.4413111,3); merge_i0j0(22.4413111,2); 
% merge_i0j0(22.44131,3); merge_i0j0(22.44131,2); 
% merge_i0j0(2.24151,2); merge_i0j0(2.24151,3); 
% merge_i0j0(2.2415101,2); merge_i0j0(2.2415101,3); 
% sfc_mode = 22.4213101; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2); merge_i0j0(sfc_mode,-2); merge_i0j0(sfc_mode,-1); merge_i0j0(sfc_mode,0);
% sfc_mode = 22.4113101; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2); merge_i0j0(sfc_mode,-2); merge_i0j0(sfc_mode,-1); merge_i0j0(sfc_mode,0);
% sfc_mode = 22.451511115; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% sfc_mode = 22.4313111; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);

% sfc_mode = 23.4414111; merge_i0j0(sfc_mode,4);

% sfc_mode = 22.4515111; merge_i0j0(sfc_mode,3);
% sfc_mode = 22.45141011; merge_i0j0(sfc_mode,2);
% sfc_mode = 22.4018111; merge_i0j0(sfc_mode,3);
% sfc_mode = 22.4718111; merge_i0j0(sfc_mode,3);
% sfc_mode = 22.4713111; merge_i0j0(sfc_mode,3);

% sfc_mode = 41.6213101; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% sfc_mode = 41.6313101; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% sfc_mode = 41.6413101; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% sfc_mode = 41.6213111; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% sfc_mode = 41.6313111; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% sfc_mode = 41.6413111; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% 

% %% asdf 
% 
% % sfc_mode = 2.2511101; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% % sfc_mode = 2.2511111; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% % sfc_mode = 2.2611101; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% % sfc_mode = 2.2611111; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2); (done i think)
% 
% % sfc_mode = 22.4513101; merge_i0j0(sfc_mode,3);
% % sfc_mode = 22.4513101; merge_i0j0(sfc_mode,2);
% % % sfc_mode = 22.4513111; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% % sfc_mode = 22.4613101; merge_i0j0(sfc_mode,3);
% % sfc_mode = 22.4613101; merge_i0j0(sfc_mode,2);
% sfc_mode = 22.4018111; merge_i0j0(sfc_mode,3); sfc_mode = 22.4018111; merge_i0j0(sfc_mode,2);
% sfc_mode = 22.4415111; merge_i0j0(sfc_mode,3);
% 
% % % % % sfc_mode = 41.6513101; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% % % % % sfc_mode = 41.6513111; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% % % % % sfc_mode = 41.6613101; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% % % % % sfc_mode = 41.6613111; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% 
% %% asdf
% 
% sfc_mode = 2.2013111; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% sfc_mode = 2.2015111; merge_i0j0(sfc_mode,3); merge_i0j0(sfc_mode,2);
% % sfc_mode = 2.2011101; merge_i0j0(sfc_mode,-1); merge_i0j0(sfc_mode,0);
% sfc_mode = 2.20141; merge_i0j0(sfc_mode,3);
% sfc_mode = 2.20141; merge_i0j0(sfc_mode,2);
% % sfc_mode = 2.2011111; merge_i0j0(sfc_mode,-1);
% % sfc_mode = 2.2011111; merge_i0j0(sfc_mode,0);
% 
% sfc_mode = 22.4015111; merge_i0j0(sfc_mode,3);
% sfc_mode = 22.4015111; merge_i0j0(sfc_mode,2);
% % sfc_mode = 22.4013101; merge_i0j0(sfc_mode,-1);
% % sfc_mode = 22.4013101; merge_i0j0(sfc_mode,0);
% sfc_mode = 22.4013111; merge_i0j0(sfc_mode,3);
% sfc_mode = 22.4013111; merge_i0j0(sfc_mode,2); %(sims incomplete)
% sfc_mode = 22.4013111; merge_i0j0(sfc_mode,-1);
% sfc_mode = 22.4013111; merge_i0j0(sfc_mode,0); %(sims incomplete)
% 
% % sfc_mode = 41.6013101; merge_i0j0(sfc_mode,3);
% % sfc_mode = 41.6013101; merge_i0j0(sfc_mode,2);
% % sfc_mode = 41.6013101; merge_i0j0(sfc_mode,-1);
% % sfc_mode = 41.6013101; merge_i0j0(sfc_mode,0);
% sfc_mode = 41.6014111; merge_i0j0(sfc_mode,3);
% sfc_mode = 41.6014111; merge_i0j0(sfc_mode,2);
% % sfc_mode = 41.6013111; merge_i0j0(sfc_mode,-1);
% % sfc_mode = 41.6013111; merge_i0j0(sfc_mode,0);
% 
% 
% % Gave error:
% % sfc_mode = 2.2013111; merge_i0j0(sfc_mode,3);
% % sfc_mode = 2.2013111; merge_i0j0(sfc_mode,2);
% % and others
% % % Error using permute_null_to_pvals>myconvert (line 106)
% % % Reference to non-existent field 'Cave1'.
% % % Error in permute_null_to_pvals>(parfor body) (line 24)
% % %         myconvert(outpath,fn,f,outpath2)
% % % Error in permute_null_to_pvals (line 23)
% % %     parfor f = 1:length(fn)
% % % Error in merge_i0j0 (line 164)
% % %         permute_null_to_pvals(sfc_mode,curr_stage)
% % % Caused by:
% % %     Reference to non-existent field 'Cave1'. 


%% Main function
function merge_i0j0(sfc_mode,curr_stage)

    %% Setup vars
    if ~exist('sfc_mode'); sfc_mode = 22.4415111; end
    if ~exist('curr_stage'); curr_stage = 2; end
    
    %% First, get list of unmerged files

    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix, ~, ~, ~, ~, ~, permutation_test] = build_sfcmode(sfc_mode, mode_subgroups);
    outpath = fullfile(getpath('path_buffer_curr'),['test3mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(curr_stage)],'i0j0_split');
    mergedpath = fullfile(outpath,'..');

    fn = dir(fullfile(outpath,'*.mat'));
    fn = {fn.name};

    %% Get list of expected final filenames (after merger)
    rawpath = fullfile(getpath('path_lfp_sample_delay'));
    fnm = dir(fullfile(rawpath,'*.mat'));
    fnm = {fnm.name};
    [~,fnmp] = cellfunu(@(f) fileparts(f),fnm);

    %% Make output folder
    %mkdir_dave(fullfile(outpath,'merged'));


    %% Search through files and keep loading until we get a sfc{1}
    for i = 1:length(fn)
        load(fullfile(outpath,fn{i}))
        if ~isempty(sfc{1})
            break
        end
    end

    %% Get list of fields for sfc{1}
    clear fieldnames1 fieldnames2 fieldnames3
    if isfield(sfc{1},'singletonfieldnames')
        fieldnames3 = sfc{1}.singletonfieldnames(:);        % If there is a field telling us the name of singleton fields, use that.
    else

                                                            % Otherwise, try to guess what they are
        sfc1 = sfc{1};
        fieldnames1 = fieldnames(sfc1);

        %% Get list of fields for sfc{~1}
        for i = 1:length(fn)
            load(fullfile(outpath,fn{i}))
            nonempty = (cellfun(@(s) ~isempty(s),sfc));
            nonempty(1)=0;
            first_nonempty = find(nonempty > 0, 1, 'first');

            if ~isempty(first_nonempty)     % We have found a struct for sfc{~1}
                break
            end

        end

        sfc2 = sfc{first_nonempty};
        fieldnames2 = fieldnames(sfc2);
        clear sfc2

        %% Get fields exclusive to sfc{1}
        j=1;
        for i = 1:length(fieldnames1)
            if sum(strcmp(fieldnames1{i},fieldnames2)) == 0
                fieldnames3{j} = fieldnames1{i};
                j=j+1;
            end
        end

    end

    %% Run the merge

    parfor i = 1:length(fnmp)
    %for i = 20
        run_the_merge(outpath,fnmp{i},fieldnames3,mergedpath)
        
    end

    %% If it's a permutation dataset, also run permute_null_to_pvals. This will shrink down the dataset size somewhat

    if permutation_test == 1
        fprintf('Detected as permutation dataset...Running permute_null_to_pvals to reduce file size \n');
        permute_null_to_pvals(sfc_mode,curr_stage)
    elseif permutation_test == 2
        fprintf('Detected as bootstrap dataset...Running bootstrap_null_to_pvals to reduce file size \n');
        bootstrap_null_to_pvals(sfc_mode,curr_stage)
        %fprintf('Detected as bootstrap dataset...merging Cavep1 and Cavep2 into dCave \n');
        %merge_distributions(sfc_mode,curr_stage);
    end
    
    %system(['rm -rf ' fullfile(outpath) ]);


end


function run_the_merge(outpath,fnmp_curr,fieldnames3,mergedpath)

        % Return without doing anything if output file exists
        outfile = fullfile(mergedpath,[fnmp_curr '.mat']);
        if exist(outfile,'file');
            fprintf('Output file already exists...skipping. \n');return;
        end

        fncurr = dir(fullfile(outpath,[fnmp_curr '*.mat']));
        fncurr = {fncurr.name};

        fprintf(['Loading ' fncurr{:} '\n'])
        s = cellfun(@(f) load(fullfile(outpath,f)),fncurr);

        clear sout
        for j = 1:length(s)
            %[i j]
            scurr = s(j).sfc;

            if ~exist('sout','var')         % If sout doesn't exist, create
                    sout = scurr;
            else                            % If it does exist, go through cells and create/append as needed

                for k = find(cellfun(@(s) ~isempty(s),scurr))

                    if length(sout) < k         % If the kth entry in sout doesn't exist, create
                        sout{k} = scurr{k};
                    else                        % If it does exist, go through fields and append as necessary

                        fld = fieldnames(scurr{k});
                        for m = 1:length(fld)
                            if ~isfield(sout{k},fld(m))                             % If the field doesn't exist, create it
                                sout{k}.(fld{m}) = scurr{k}.(fld{m});
                            elseif (sum(strcmp(fld{m},fieldnames3)) == 0 || strcmp(fld{m},'adj_was_missing'))             % If the field is a sfc{1} field, skip it; otherwise, append.
                                sout{k}.(fld{m}) = cat(2,sout{k}.(fld{m}),scurr{k}.(fld{m}));
                            end
                        end
                    end
                end
            end
        end
        if ~isempty(s)
            sfc = sout;
            save(outfile,'sfc','-v7.3');
        end
end


