
function os = get_options(sfc_mode, curr_stage,i0,j0,outname_embed_i0j0,file_list,filenum)
    %% function options_struct = get_options(sfc_mode,curr_stage)

    if nargin < 2
        error('Need at least 2 inputs - sfc_mode and curr_stage.')
    end
    
    if ~exist('file_list','var'); file_list = []; end
    if ~exist('filenum','var'); filenum = []; end
    
    if ~isempty(file_list) && ~isempty(filenum)
        fname = file_list{filenum};
    end
    
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode, do_partial_coherence, ctgsetli_mode2, Nwind_mode] = build_sfcmode(sfc_mode, mode_subgroups);
    sfc_mode_group = floor(sfc_mode); % Mode group might be more than 1 dimension, so use floor instead (i.e. sfc_mode = 23.2 would yield mode_group == [2 3] instad of 23)
    
    skip_Jackknife_for_ctg_alltrials17 = 1;

    exclude_clipping_trials=1;
    
    if sfc_mode_group == 6     % Raster plotting mode
        curr_stage = 5; % For looking at the whole set of activity.
    end
    
    if sfc_mode_group == 3 || sfc_mode_group == 23 || sfc_mode_group == 33 || sfc_mode_group == 45
        is_spectrogram = 1;
    else
        is_spectrogram = 0;
    end
    
    dt = get_dt;
    
    
    
    thresh = 5;                 % Firing rate threshold for rejecting traces
    switch tapers_mode
        case 0; params.tapers = get_tapers;
        case 3; params.tapers = [1,1];
        case 4; params.tapers = [3,5];
        case 5; params.tapers = [5,9];
        case 7; params.tapers = [10,19];
        case 8; params.tapers = [2,3];
        case 9; params.tapers = [20,39];
        otherwise; params.tapers = get_tapers;
    end
    
    params.fpass = [0 120];
    params.err = [1 0.05];      % Using asymptotic method.
    
    if tapers_mode >= 1 && tapers_mode <= 2
        tapers_T = (diff(get_stagesir(curr_stage))+1)/1e3;
        tapers_desired_bandwidth = 20;
        tapers_W = tapers_desired_bandwidth/2;      % 2W is the bandwidth (minf to maxf; verified in Chronux using sinusoidal test signal) 
        tapers_TW = round(tapers_T * tapers_W );
        params.tapers = [tapers_TW 2*tapers_TW-1];
    end
    if tapers_mode == 2                 % Jackknife
        params.err(1) = 2;
    end
    
    params.Fs = 1/dt;
    params.trialave = 1;
    
    tapers_mode_orig = tapers_mode;     % These are used later, only for the skip_jackknife switch; it's fine, really.
    params_err_orig = params.err(1);
    
    needed_spikes = get_needed_spikes(ue_pairs);
    
    needed_lfp = get_needed_lfp(ue_pairs);
    
    
    
    % Set up Nwind & estimate bandwidths
    if is_spectrogram
        switch Nwind_mode
            case 0
                Nwind = round(0.3/dt);
            case 1
                Nwind = round(0.6/dt);
            case 2
                Nwind = round(0.8/dt);
            case 3
                Nwind = round(1.0/dt);
        end
        
        f_fullbandwidth = 2 * params.tapers(1) / (Nwind*dt); % Recall tapers(1) = TW, therefore W = tapers(1)/T. Thus, 2W=2*tapers(1)/T is min to max bandwidth
        
        
    elseif sfc_mode_group == 52 % Set sgolay window length for unit analysis
        
            switch tapers_mode
                case 0; Nwind= diff(get_stagesir(curr_stage))+1; Nwind = Nwind * dt;    % No smoothing; just average
                case 1; Nwind = 0.151;
                case 2; Nwind = 0.301;    % This is the full width - NOT THE HALF WIDTH
                case 3; Nwind = 0.601;
                case 4; Nwind = 0.751;
                case 5; Nwind = 0.051;
            end
            Nwind = round(Nwind/dt); % in samples
            
            f_fullbandwidth = Nwind*dt;
            
    else                        % Estimate window size for FFC spectra.
        Nwind= diff(get_stagesir(curr_stage))+1;
        f_fullbandwidth = 2 * params.tapers(1) / (Nwind*dt); % Recall tapers(1) = TW, therefore W = tapers(1)/T. Thus, 2W=2*tapers(1)/T is min to max bandwidth
    end
        
    fract_overlap = 0.9;
    
    os = vars_push;
    
    if nargin > 2
        [outname,outpath] = get_outname (os);
        os.outname = outname;
        os.outpath = outpath;
    end


end




%% Helper functions
function out = get_needed_spikes(ue_pairs)
    % Returns true if spike traces need to be calculated
    out = any(ue_pairs == [0,2,3,5,7]);
end

function out = get_needed_lfp(ue_pairs)
    % Returns true if lfp traces need to be calculated
    out = any(ue_pairs == [0,2,3,4,6]);
end


function [outname,outpath] = get_outname (os)
    sfc_mode = os.sfc_mode;
    curr_stage = os.curr_stage;
    outname_embed_i0j0 = os.outname_embed_i0j0;
    filename_orig = os.fname;
    fname_suffix = os.fname_suffix;

    % Generate outpath
    fprintf(['Chosen SFC mode: ' fname_suffix ' \n']);
    outpath = fullfile(getpath('path_buffer_curr'),['test2mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(curr_stage)]);
    if outname_embed_i0j0
        outpath = fullfile(outpath,'i0j0_split'); [~,filename,~] = fileparts(filename_orig);
        filename = [filename '_i' num2str(i0(1),'%2.2d') '_j' num2str(j0(1),'%2.3d') '.mat'];
    else
        filename = filename_orig;
    end
    %mkdir_dave(outpath);
    outname = fullfile(outpath,filename);
end