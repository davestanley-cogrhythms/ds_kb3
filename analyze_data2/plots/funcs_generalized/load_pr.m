
function [wrkspc_buffer, out] = load_pr(wrkspc_buffer,sfc_mode,curr_stage_sp,freqband_perm,bad_any,opts_perm,opts_exclude)


%% Setup basic parameters
% clc
addpath(genpath('../funcs_plots_preferred/'));

% Data loading parameters
recalc_permute = 0;
reload_data = 0;
if get_iscromer; stage_range = [2,3];
else stage_range = [2,3];
end

    % SFC mode
    if ~exist('sfc_mode','var'); sfc_mode = 22.4013111; end

    % Stage selection
    if ~exist('curr_stage_sp','var'); curr_stage_sp = 2; end
    
    % Set up other parameters
    if ~exist('freqband_perm','var'); freqband_perm = [16 20]; end
    if ~exist('bad_any','var');
        bad_any = [];           % Assuming nothing bad, unless otherwise specified
    end
    
        
    
    
    % Test plot switches
    plot_sums = 0;
    plot_imagesc_passes = 0;
    

if sfc_mode == 5
    % If sfc_mode is 5, load preferred units (based on firing rates) and
    % skip the stuff below.
    
    do_bh = 1;
    [sp,sr,ss,pr_fnames,pr_funames1D] = extract_spsr(curr_stage_sp, bad_any,do_bh);
    
    % Package output
    out.sig_cells = sp;
    out.Cave1 = [sr(:,1:2:end)];
    out.Cave2 = [sr(:,2:2:end)];
    out.mypairs = [[1:size(sp,1)]' [1:size(sp,1)]'];


else
    % If sfc_mode is not 5, load permutation tests as normal
    
    %% Derive some basic things based on parameters

    % Derive suffix name and switches
    [~, sfc_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix_sfc, do_adjacent, ctgsetli_mode, thinning_mode, tapers_mode, baseline_subtract, permutation_test, ue_pairs, coh_debias_mode] = build_sfcmode(sfc_mode, sfc_subgroups);
    fprintf(['Chosen SFC mode: ' fname_suffix_sfc ' \n']);
    
    
    os = get_options(sfc_mode, curr_stage_sp);



    %% Load data    

    % Estimate file list and file_ranges
    [~, file_range] = estimate_filelist(sfc_mode, curr_stage_sp);

    path_matfiles_store = getpath('path_matfiles_store');


    % Create buffer variable
    if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end

    % Characterize data type;
    data_type = mode2datatype(sfc_mode);
    is_spectrogram = get_is_spectrogram(sfc_mode);


    % Load the Perm data if it is in perm mode
    perm_mode = sfc_mode;
    buff_mode = perm_mode;
    buff_fieldname = ['per_' mode2modename(buff_mode)];
    buff_func = @() func_recalc_perms(buff_mode,stage_range,file_range);
    wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_permute,reload_data);



    % Pull in the data
    vars_pull(wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sp)));
    % if strcmp(data_type,'FFC');
    vars_pull(wrkspc_buffer.(buff_fieldname).pairsdat.(stagename(curr_stage_sp)));
    % end


    mypairs = mypairs(:,:,end)';
    Ncells_shift = Ncells_shift(:);
    Ncells_shift = Ncells_shift(:);
    Nelects_shift = Nelects_shift(:);
    lumap=lumap(:);
    
    % Set up mypairs
    if strcmp(data_type,'FFC') || strcmp(data_type,'PSD')
        mypairs = mypairs + [Nelects_shift Nelects_shift];
    elseif strcmp(data_type,'SFC') || strcmp(data_type,'FR')
        mypairs = mypairs + [Ncells_shift Ncells_shift];
    end
    
    % For FR time series, shift f so t=0 is aligned with sample on
    if strcmp(data_type,'FR')
        temp = get_stagesir(curr_stage_sp);
        f = f + temp(1)/1e3;
    end

    % For spectrogram, center f2 time series vector
    if is_spectrogram
        %temp = get_stagesir(curr_stage_sfc);
        %f2 = f2 + temp(1)/1e3;
        if exist('stagesir','var')
            if any(f2 > 100)                    % If any time series are greater than 100, we know its in milliseconds.
                f2 = f2/1e3 + stagesir(1)/1e3;  % This works with most recently calculated spectrograms using precalc
            else
                f2 = f2 + stagesir(1)/1e3;      % Annoying case for old spectrogram data
            end
        else
            f2 = f2 - 2.0; warning('Stagesir not found; guessing time range of spectrogram.');
        end
    else
        f2 = [];
    end
    
    % Load bads specific to pr
    if exist('opts_exclude','var')
        [bad_any_sp] = load_bads(sfc_mode,curr_stage_sp,opts_exclude,wrkspc_buffer.currmd.md,mypairs,funames1D);
    else
        error('Update to code - must include opts_exclude. Change load_pr(wrkspc_buffer,sfc_mode,curr_stage_sp,freqband_perm,bad_any,opts_perm) to load_pr(wrkspc_buffer,sfc_mode,curr_stage_sp,freqband_perm,bad_any,opts_perm,opts_exclude)');
        bad_any_sp = bad_any;
    end

    %% Pull out the permutation test data

    s = wrkspc_buffer.(buff_fieldname).(stagename(curr_stage_sp));
    
    s.f = f;
    s.f2 = f2;
    
    [sig_cells] = plot_permute_scatter_fast(s, bad_any_sp,freqband_perm, opts_perm,is_spectrogram);



    %% Run_diagnostic
    if ~isempty(bad_any_sp)
        if length(bad_any_sp) == size(sig_cells,1)     % They must be the same size to run the test. If not, skip.
            %% Bad_any test
            % Make sure that all the bad_any's are actually being excluded!
            sum_goods = sum(sig_cells);
            sum_goods_and_not_bad = sum(sig_cells & repmat(~bad_any_sp(:),1,size(sig_cells,2)));
            if any(sum_goods ~= sum_goods_and_not_bad);
                warning('Bads are not being properly removed!');
            end
        end
    end

    %% Plotting
    if plot_sums
        %% Plot_sums
        figure; bar(sum(sig_cells));

    end

    if plot_imagesc_passes
        %% Plot_imagesc_passes
        figure; imagesc(sig_cells);

    end

    %% Package output
        out.sig_cells = sig_cells;
        out.Cave1 = Cave1;
        out.Cave2 = Cave2;
        out.mypairs = mypairs;
        out.funames1D = funames1D;
        out.bad_any_sp = bad_any_sp;
        out.abscissa = s.f;
        out.os = os;

end

end
