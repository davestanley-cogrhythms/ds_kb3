

function run_batch
% Matlab version of run batch script. Makes for much easier coding than
% bash!

    addpath(fullfile('plots','funcs_plots_preferred'));

    %run_allstages(41.6013101);
%     run_allstages(2.2005111);
%     run_allstages(2.201511101);
%     run_allstages(2.201501101);
%     run_allstages(2.001511101);
%     run_allstages(2.001501101);
%     run_allstages(2.301711101);
%     run_allstages(2.371711101);

%     run_allstages(22.441511102);

%     run_allstages(22.461810303);
%     run_allstages(22.461310303);
%     run_allstages(22.461810302);
%     run_allstages(22.461310302);

%     run_allstages(23.4513101031);
%     run_allstages(23.4513101033);
    
%     run_allstages(22.471810102);
%     run_allstages(22.471810103);
%     run_allstages(22.481810102);
%     run_allstages(22.481810103);
    
%     run_allstages(23.471810103);
    
    
%     
%     run_allstages(45.6013111012);
    run_allstages(45.6013113041);
    run_allstages(23.4013113041);
%     run_allstages(45.6013103041);
%     run_allstages(23.4013103041);
%     run_allstages(23.4013111013);
%     run_allstages(23.4018111013);
    
    

    % run_allstages(22.471711111);
%     run_allstages(22.471511101);
    
%     run_allstages(22.491510310);
%     run_allstages(22.491710310);
%     run_allstages(22.491510311);
%     run_allstages(22.491710311);
%     
%     run_allstages(22.451511103);
%     run_allstages(22.451711103);
%     run_allstages(22.401711101);


%     run_allstages(22.401511101);
    
%     run_allstages(52.740001002);
%     run_allstages(52.750001003);
    
    run_allstages(52.700200004);
%     run_allstages(52.750201003);

    
%     run_allstages(2.241711102);
%     run_allstages(2.251711103);
%     run_allstages(2.241711100);
%     run_allstages(2.251711100);
%     run_allstages(3.201711101);
%     run_allstages(3.201911101);
%     run_allstages(23.4013111);
%     run_allstages(23.4014111);
%     run_allstages(52.7402010);
%     run_allstages(52.7403010);
%     run_allstages(52.7502010);
%     run_allstages(52.7503010);
%       run_allstages(22.401311101);

%       run_allstages(22.441311102);
%       run_allstages(22.441411102);
%       run_allstages(22.441511102);
%       run_allstages(22.441311103);
%       run_allstages(22.441411103);
%       run_allstages(22.441511103);
% 
%       run_allstages(22.451311102);
%       run_allstages(22.451411102);
%       run_allstages(22.451511102);
%       run_allstages(22.451311103);
%       run_allstages(22.451411103);
%       run_allstages(22.451511103);
%     run_allstages(22.451411103);
%       


%       run_allstages(23.441410102);
%       run_allstages(23.441410103);
%       run_allstages(23.451410102);
%       run_allstages(23.451410103);
      
%       run_allstages(22.491410101);
%       run_allstages(52.790200001);
%       run_allstages(52.700201001);
%     run_allstages(23.491410301);
      %run_allstages(22.4515111);
%     run_allstages(22.4718111);
%     run_allstages(22.4713101);
%     run_allstages(22.4517101);
%     run_allstages(2.0015111);
%     run_allstages(22.4918103);
%     run_allstages(52.7001000);
%     run_allstages(52.7400010);
%     run_allstages(52.7501000);
%     run_allstages(52.7900000);
%     run_allstages(52.7005010);
%     run_allstages(52.770001000);
%     run_allstages(52.770001001);
%     run_allstages(52.770201000);
%     run_allstages(52.770201001);
%     run_allstages(52.7401010);

%     run_allstages(2.271511101);
%     run_allstages(2.271911101);
    %run_allstages(2.2415100);
    
    %run_allstages(22.40171111);
    %run_allstages(22.45171111);
    %run_allstages(22.40141011);
%     run_allstages(41.6013111);
%     run_allstages(41.6515111);
%     run_allstages(41.6015111);
    %run_allstages(22.4514111); 
    %run_allstages(22.4014111);

end

function run_allstages (sfc_mode)

%     total_stages = 2;
%     run_qsub(sfc_mode,2,total_stages);
%     run_qsub(sfc_mode,3,total_stages);

    
    total_stages = 1;
    run_qsub(sfc_mode,4,total_stages);
    

end


function run_qsub(sfc_mode,curr_stage,total_stages)
    % Parameters
    server_mode = 2; % Server mode:
                    % 1-Run on server using compiled code;
                    % 2-Run on server using Matlab code
                    % 3-Run locally using Matlab code
                    % 4-Run using GPU's
    i0j0_split = 0;
        % 0 - Do 1 whole file per simulation (this saves time on loading
        %     lfp/spike traces.
        % 1 - Do a total of "desired_number_of_nodes" simulations, with
        %     pairs evenly distributed throughout.
    
    desired_number_of_nodes = 400; % (Only for sim_mode 2)
    desired_number_of_nodes = round(desired_number_of_nodes/total_stages);
    
    % Decode sfc_mode
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2] = build_sfcmode(sfc_mode, mode_subgroups);
    sfc_mode_group = floor(sfc_mode); 
    mode_firstdecimal = mode_subgroups(1);
    mode_sixthdecimal = mode_subgroups(6);
    
    do_precalc = 1;
    if do_partial_coherence; do_precalc = 0; end
    
    
    % Load Metadata struct
    path_matfiles_store = getpath('path_matfiles_store');
    wrkspc_buffer = struct;

    buff_fieldname = 'currmd';
    buff_func = @() func_recalc_md();
    wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,0,0);
    md = wrkspc_buffer.(buff_fieldname).md;
    md=md(1);

    [ctgsetli] = get_ctgsetli(sfc_mode,md,ctgsetli_mode,ctgsetli_mode2,1);

    % Estimate i0 range
    i0_unpaired_str=['1:' num2str(size(ctgsetli,2))];
    i0_paired_str = ['1:' num2str(floor(size(ctgsetli,2)/2)*2)];
    %i0_paired_str = [num2str(round(size(ctgsetli,2)/2)+1) ':' num2str(floor(size(ctgsetli,2)/2)*2)];    % Correction
    %i0_paired_str = [num2str(round(size(ctgsetli,2)/2)+1) ':' num2str(floor(size(ctgsetli,2)/2)*2)];    % Correction
%     i0_unpaired_str=['1:10'];
%     i0_paired_str = i0_unpaired_str;

    % Estimate desired runtime
    if ~do_partial_coherence
        if ~permutation_test; myhours = '1';
        else myhours = '36';
        end
        if server_mode == 4
            if ~permutation_test; myhours = '1';
            else myhours = '5';
            end
        end
    else
        if ~permutation_test; myhours = '3';
        else myhours = '196';
        end
    end
    
    
    
    if i0j0_split == 0
        % Estimate list of files
        path_md = fullfile(getpath('path_metadata'),'sansunits');
        file_list = get_filelist(path_md);
    elseif i0j0_split == 1
        % Get full list of simulations
        load Nall.mat
        switch mode_firstdecimal
            case 2; Npairs_arr = Nunits;
            case 3; Npairs_arr = Nunits_field;
            case 4; Npairs_arr = Nfield_field;
            case 5; Npairs_arr = Nunits_units;
            case 6; Npairs_arr = Nlfp;
        end

        % Generate array of simulation indicies and corresponding file numbers
        % and pair numbers
        pairN = arrayfunu(@(x) 1:x, Npairs_arr);
        fileN = arrayfunu(@(x,y) x*ones(1,y), 1:length(Npairs_arr),Npairs_arr);
        pairN = horzcat(pairN{:});
        fileN = horzcat(fileN{:});
        Nsimulations = length(pairN);
        sims_per_node = ceil(Nsimulations / desired_number_of_nodes);
        actual_number_of_nodes = ceil(Nsimulations / sims_per_node);
    end

    
    % Setup i0 depending on whether we're in permute mode or not
    if mode_sixthdecimal == 0; i0_str = i0_unpaired_str;
    else i0_str = i0_paired_str;
    end
    
    switch server_mode
        case 1
            if i0j0_split == 0
                for ind = 1:length(file_list)
                    fprintf(['Running file ' num2str(ind) '\n']);
                    error('Verify this works with precalc');
                    system(['qsub -l h_rt=' num2str(myhours) ':30:00 run_standalone_job.sh run_analysis_cells_binwrapp ' num2str(sfc_mode,'%10.12f') ' ' num2str(curr_stage) ' ' num2str(ind) ' ' num2str(ind) ' "' i0_str '" ' num2str(i0j0_split) ' ' num2str(do_precalc) '']);
                    fprintf('\n');
                end
            elseif i0j0_split == 1
                for ind = 1:actual_number_of_nodes
                    pairsstart = (ind-1)*sims_per_node + 1;
                    pairsstop = min(Nsimulations,ind*sims_per_node);
                    fprintf(['Running pairs ' num2str(pairsstart) ' to ' num2str(pairsstop) '\n']);
                    error('Verify this works with precalc');
                    system(['qsub -l h_rt=' num2str(myhours) ':30:00 run_standalone_job.sh run_analysis_cells_binwrapp ' num2str(sfc_mode,'%10.12f') ' ' num2str(curr_stage) ' ' num2str(pairsstart) ' ' num2str(pairsstop) ' "' i0_str '" ' num2str(i0j0_split) ' ' num2str(do_precalc) '']);
                    fprintf('\n');
                end
            end
            
        case 2
            if i0j0_split == 0
                for ind = 1:length(file_list)
                    fprintf(['Running file ' num2str(ind) '\n']);
                    system(['qsub -l h_rt=' num2str(myhours) ':30:00 -l mem_total=95G matlab_single_node_batch.sh "setup_paths_n_run(@run_analysis_cells,' num2str(sfc_mode,'%10.12f') ',' num2str(curr_stage) ',' num2str(ind) ',' num2str(ind) ',' i0_str ',' num2str(i0j0_split) ',' num2str(do_precalc) ')" localOutput']);
                    fprintf('\n');
                end
            elseif i0j0_split == 1
                for ind = 1:actual_number_of_nodes
                    pairsstart = (ind-1)*sims_per_node + 1;
                    pairsstop = min(Nsimulations,ind*sims_per_node);
                    fprintf(['Running pairs ' num2str(pairsstart) ' to ' num2str(pairsstop) '\n']);
                    system(['qsub -l h_rt=' num2str(myhours) ':30:00 matlab_single_node_batch.sh "setup_paths_n_run(@run_analysis_cells,' num2str(sfc_mode,'%10.12f') ',' num2str(curr_stage) ',' num2str(pairsstart) ',' num2str(pairsstop) ',' i0_str ',' num2str(i0j0_split) ',' num2str(do_precalc) ')" localOutput']);
                    fprintf('\n');
                end
            end
        case 3
            if i0j0_split == 0
                for ind = 1:length(file_list)
                    fprintf(['Running file ' num2str(ind) '\n']);
                    setup_paths_n_run(@run_analysis_cells,sfc_mode,curr_stage,ind,ind,i0_str,i0j0_split, do_precalc);
                end
            elseif i0j0_split == 1
                for ind = 1:actual_number_of_nodes
                    pairsstart = (ind-1)*sims_per_node + 1;
                    pairsstop = min(Nsimulations,ind*sims_per_node);
                    fprintf(['Running pairs ' num2str(pairsstart) ' to ' num2str(pairsstop) '\n']);
                    setup_paths_n_run(@run_analysis_cells,sfc_mode,curr_stage,pairsstart,pairsstop,i0_str,i0j0_split, do_precalc);
                end
            end
        case 4
            for ind = 1:length(file_list)
                fprintf(['Running file ' num2str(ind) '\n']);
                system(['qsub -l h_rt=' num2str(myhours) ':30:00 -l gpus=1 matlab_single_node_batch.sh "setup_paths_n_run(@run_analysis_cells,' num2str(sfc_mode,'%10.12f') ',' num2str(curr_stage) ',' num2str(ind) ',' num2str(ind) ',' i0_str ',' num2str(i0j0_split) ',' num2str(do_precalc) ')" localOutput']);
                fprintf('\n');
            end
    end
    
end

