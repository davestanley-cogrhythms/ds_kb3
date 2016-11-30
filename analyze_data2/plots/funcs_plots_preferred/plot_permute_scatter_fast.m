
function [sig_cells] = plot_permute_scatter_fast(s, bad_any,freqband_perm, opts,is_spectrogram)
    %% plot_SchAvsB_scatter - Plot significant differences AMPLITUDE across SchA and SchB
    
    if ~exist('opts','var'); opts = Opts_Perm; end
    if ~exist('is_spectrogram','var'); is_spectrogram = 0; end
    
    % Pull in options
    do_bh0 = opts.do_bh0;
    do_phi = opts.do_phi;
    split_plusminus = opts.split_plusminus;
    alpha0 = opts.alpha0;
    alpha_bh0 = opts.alpha_bh0;
    do_quantiles_mode = opts.do_quantiles_mode;
            chosen_quantile = opts.chosen_quantile;
            upper_quantile = opts.upper_quantile;
    timeband_perm = opts.timeband_perm;
    
    % Setup defaults
    %if isempty(sens_chosen_condition); sens_chosen_condition = 1:size(s.pvals_Cave,3); end
    
    chosen_cells = ~bad_any;
    
    alpha = alpha0;
    if do_bh0; alpha = alpha_bh0; end
    
    % Get permmode 1 data for adding to graph below

    if ~do_phi; pvals_Cave_perm = s.pvals_Cave(:,:,:);
    else pvals_Cave_perm = s.pvals_phi(:,:,:); end
    f = s.f;
    f2 = s.f2;
    
    if isempty(freqband_perm); error('Freqband_perm not assigned.'); end
    if is_spectrogram && isempty(timeband_perm); error('Timeband_perm not assigned.'); end
    
    if ~do_quantiles_mode
        if is_spectrogram; pvals_Cave_perm = unwrap_3Dpls(pvals_Cave_perm,f,f2); end
        [sig_cells, pvals1, p0] = calc_significant_permutation2D(f,f2,pvals_Cave_perm,alpha,bad_any,mean(freqband_perm),do_bh0,is_spectrogram,mean(timeband_perm));
    else
%         [sig_cells] = calc_significant_quantiles(f,pvals_Cave_perm,alpha,bad_any,mean(freqband_perm),do_bh0);

        if ~do_phi; dCave = s.Cave1 - s.Cave2;
        else dCave = abs(s.phi1 - s.phi2);
        end
        if is_spectrogram; dCave = unwrap_3Dpls(dCave,f,f2); end
        [sig_cells] = calc_significant_quantiles(f,f2,dCave,chosen_quantile,upper_quantile,bad_any,mean(freqband_perm),is_spectrogram,mean(timeband_perm));
        
    end
%     
    %% Manipulate sig_cells as needed
    % Split into positive and negative directions
    
    switch split_plusminus
        case 0
            % do nothing
        case 1
            % Include both positive and negative
            sig_cells = do_split_plusminus(sig_cells,s,mean(freqband_perm),do_phi,is_spectrogram,mean(timeband_perm));
        case 2
            % Include only positive
            sig_cells_temp = do_split_plusminus(sig_cells,s,mean(freqband_perm),do_phi,is_spectrogram,mean(timeband_perm));
            sig_cells = sig_cells_temp(:,1:2:end);   % Take odds
            clear sig_cells_temp
        case 3
            % Include only negative
            sig_cells_temp = do_split_plusminus(sig_cells,s,mean(freqband_perm),do_phi,is_spectrogram,mean(timeband_perm));
            sig_cells = sig_cells_temp(:,2:2:end);  % Take evens
            clear sig_cells_temp
    end
    
end



function [sig_cells, pvalsf, p0] = calc_significant_permutation2D(f,f2,pvals,alpha,bad_any,freqval,do_bh_stepup,is_spectrogram,timeband_perm)
    %% function calc_significant_permutation2D
    % Newer version of calc_significant_permutation, which works with 2D
    % data. Not well tested.

    if ~exist('do_bh_stepup','var'); do_bh_stepup = 1; end      % Use BH setup up algorithm to enforce fixed false discovery rate, rather than just simply doing repeated t test
    plot_on = 0;
    plot_debug = 0;
    
    if length(bad_any) ~= size(pvals,2)
        warning('Size of bad_any doesnt match number of elements in pvals. Disregarding bad_any');
        bad_any = [];
    end
    
    if isempty(bad_any); bad_any = false(size(pvals,2),1); end
    
    
    % Look only at a specific frequency
    if ~is_spectrogram
        index = find_closest(f,freqval);
        pvalsf = squeeze(pvals(index,:,:));
    else
        index = find_closest(f,freqval);
        index2 = find_closest(f2,timeband_perm);
        pvalsf = squeeze(pvals(index,:,:,index2));
    end
    
    
    if do_bh_stepup
        %[sig_cells, p0] = bh_step1D(pvalsf,alpha,bad_any);    % My code
        [sig_cells, p0] = bh_step(pvalsf,alpha,bad_any);    % My code
        %sig_cells2 = fdr_bh(pvals,alpha);  % Mathworks code - they seem to give the same result
    else
        sig_cells = pvalsf < alpha;
        p0 = alpha;
    end
    
    bad_any2 = repmat(bad_any(:),1,size(sig_cells,2));
    sig_cells = sig_cells & ~bad_any2;
    
    
end

function sig_out = do_split_plusminus(sig_cells,s,freqband_perm,do_phi,is_spectrogram,timeband_perm)
    %% function do_split_plusminus
    % Calculate difference statistic
    if ~do_phi
        dCave = s.Cave1-s.Cave2;
    else
        dCave = s.phi1-s.phi2;
    end
    
    
    if ~is_spectrogram
        dCave_stats = calc_pls_stats(s.f,dCave,freqband_perm,'do_mean_ctgs',0);
    else
        dCave = unwrap_3Dpls(dCave,s.f,s.f2);
        index = find_closest(s.f,freqband_perm);
        index2 = find_closest(s.f2,timeband_perm);
        dCave_stats = squeeze(dCave(index,:,:,index2));
    end
    
    if do_phi; dCave_stats = get_circ_diff(dCave_stats); end
    
    

    sz = size(sig_cells);
    sig_out = false(sz(1),sz(2)*2);
    for i = 1:sz(2)
        sig_out(:,i*2-1) = sig_cells(:,i) & dCave_stats(:,i) >= 0;  % Entries 1,3,5,... - Significant + diff is positive
        sig_out(:,i*2) = sig_cells(:,i) & dCave_stats(:,i) < 0;  % Entries 1,3,5,... - Significant + diff is negative
    end
end


function [sig_cells] = calc_significant_quantiles(f,f2,dCave,chosen_quantile,upper_quantile,bad_any,freqval,is_spectrogram,timeband_perm)
    %% [sig_cells] = calc_significant_quantiles(f,dCave,chosen_quantile,upper_quantile,bad_any,freqval)
    quantile_mode = 1;
            % Mode 1 - evaluate at quantiles of each ctg individually
            % Mode 2 - pool all ctgs together and then evaluate quantiles
            
    switch quantile_mode
        case 1
            Q = quantile(dCave(:,~bad_any,:,:),chosen_quantile,2); Q = repmat(Q,[1,size(dCave,2),1,1]);
        case 2
            temp = dCave(:,~bad_any,:,:,:);
            temp = permute(temp,[1,4,2,3]);       % Swap data so it's: (freq,time,cells,ctgs)
            temp = temp(:,:,:);                   % Collapse cells and ctg data into dimension 3
            Q = quantile(temp,chosen_quantile,3); % Take quantil
            Q = repmat(Q,[1,1,size(dCave,2),size(dCave,3)]);    % Re-expand electrode (dim 2) and ctg data (dim 4)
            Q = permute(Q,[1,3,4,2]);
    end
    
    if upper_quantile
        sig_cells = dCave > Q;
    else
        sig_cells = dCave < Q;
    end
    
%     index = find_closest(f,freqval);
%     sig_cells = squeeze(sig_cells(index,:,:,:,:));
    
    % Look only at a specific frequency
    if ~is_spectrogram
        index = find_closest(f,freqval);
        sig_cells = squeeze(sig_cells(index,:,:,:,:));
    else
        index = find_closest(f,freqval);
        index2 = find_closest(f2,timeband_perm);
        sig_cells = squeeze(sig_cells(index,:,:,index2));
    end
    
    
end
