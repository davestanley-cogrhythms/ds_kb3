
function [sig_cells, p0] = bh_step1D(pvals0,qval, bad_any)
    % Old, depreceated function that only works with 1D

    plot_on = 0;

    if ~exist('bad_any','var'); bad_any = []; end
    
    if ~isempty(bad_any)        % If we're excluding bad cells, only consider good cells for estimation of p0
        pvals = pvals0(~bad_any);
    else
        pvals = pvals0;
    end

    if isvector(pvals); pvals = pvals(:); end
    
    if ~isvector(pvals0); warning('Using 2D data for bh_step that only works with 1D'); end
    
    sig_cells = false(size(pvals));
    p0 = zeros(1,size(pvals,2));
    
    if plot_on; figure; end
    
    for i = 1:size(pvals,2)

        pvals_curr = pvals(:,i);
        pvalss = sort(pvals_curr);
        m = length(pvalss);
        ind = find( pvalss ./ ((1:m)'*qval/m) <= 1, 1, 'last');
        p0_curr = pvalss(ind);

        if ~isempty(p0_curr)
            sig_cells(:,i) = pvals_curr <= p0_curr;
        else
            sig_cells(:,i) = false(size(pvals_curr));
            p0_curr = 0;
        end
        
        p0(i) = p0_curr;
        
        if plot_on
            clf;
            plot(1:m,pvalss);
            hold on; plot(1:m,((1:m)'*qval/m),'r');
            plot(ind,p0_curr,'k*');
            legend('pValues','Threshold','p0');
            
        end
    end
    
    if ~isempty(bad_any)            % Then use this estimation of p0 to pick out significant cells.
        sig_cells = pvals0 <= p0;
    end

end