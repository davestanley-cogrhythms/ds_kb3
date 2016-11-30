

function [sig_cells, pvalsf, p0] = calc_significant_permutation(f,pvals,alpha,bad_any,freqval,do_bh_stepup)

    if ~exist('do_bh_stepup','var'); do_bh_stepup = 1; end      % Use BH setup up algorithm to enforce fixed false discovery rate, rather than just simply doing repeated t test
    plot_on = 0;
    plot_debug = 0;
    
    
    % Look only at a specific frequency
    %index = find(f >= freqval, 1, 'first');
    index = find_closest(f,freqval);
    pvalsf = squeeze(pvals(index,:));
    
    if do_bh_stepup
        
        [sig_cells, p0] = bh_step1D(pvalsf,alpha,bad_any);    % My code
        %sig_cells2 = fdr_bh(pvals,alpha);  % Mathworks code - they seem to give the same result
        
        
    else
        sig_cells = pvalsf < alpha;
        p0 = alpha;
    end
    
    
    
end


%%