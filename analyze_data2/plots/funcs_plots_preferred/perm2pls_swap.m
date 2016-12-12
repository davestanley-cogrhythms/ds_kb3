function [pls_new] = perm2pls_swap(s_perm,perm2pls_dophi,perm2pls_do_bh,bad_any,return_mode,allow_signed,split_plusminus)
    alpha = 0.05;
    
    if ~exist('allow_signed','var');    % If 0, returned pls is always positive; otherwise it can be +ve or -ve
        allow_signed = 0;               % depending on sign of original Cave1-Cave2
    end
    
    if ~exist('split_plusminus','var');    
        split_plusminus = 0;
                    % 1-Return everything; 1-return pos+negative inseparate columns;
                    % 2-return only positive (others set to 0);
                    % 3-return only negative;
                    % 4-return only significant (only makes sense for return_mode=4)
    end
    
    if ~exist('return_mode','var'); 
        return_mode = 4;
            % 1 returns pass/fail values of significance test for each cell
            % 2 returns area under cdf for each cell (from center of distribution; dFFC=0 corresponds to pls=0; dFFC=Inf corresponds to pls=1)
            % 3 returns z-score metric, similar to Buschmann
            % 4 uses actual z-score (should be saved directly by
            %               permute_null_to_pvals in later versions of code).
    end
        
    if return_mode > 1; perm2pls_do_bh=0; end
    
    if ~perm2pls_dophi; pvals = s_perm.pvals_Cave;
    else pvals = s_perm.pvals_phi;
    end
    
    % Summate cells passing test
    % ~~~~~~~~~~~~~~ THIS NO LONGER WORKS, NEED PLS TO BE SIZE NCELLS ~~~~~~~~~~~`
    %pls_new = sum(sig,2)/sum(~bad_any) *100;     % For percent electrodes
    %pls_new = sum(sig,2);                        % For total
    
    % Calculate number of significant cells if needed
    switch return_mode
        case 1
            % Get significant cells
            sig = get_sig(pvals,perm2pls_do_bh,alpha,bad_any);
            
            %pls_new = sig*100;                              % For percent of cells
            pls_new = sig*sum(~bad_any);                    % For number of cells
            
            % Do split_plusminus as needed
            renormalize_on = 0;
            pls_new = do_splitplusminus(pls_new,sig,split_plusminus,s_perm,perm2pls_dophi,renormalize_on);

  
        case 2
            pls_new = 1-pvals;                    % For area under CDF starting from center (based on pvalue)
        case 3
            pls_new = norminv(1-pvals/2,0,1);   % This transforms the pvals to dFFC, assuming dFFC ~ N(0,1) (i.e. normal distribution).
            pls_new=fix_infs(pls_new);          % This fixes the infs issue. It arises because we still have some pvals that are equal to 0; this is an unrealistic side effect
                                                % of the permutation test and produces infs in the norminv transform. The function fix_infs
                                                % converts these 0 pvals to estimates of what they could be by enforcing that dFFC be normally
                                                % distributed (seems to be reasonable since dFFC is close to normal with the infs removed).
            % Pls_new is now abs(N(0,1)) (since all values are positive).
            % We want to convert this to a z-score. Therefore normalize to
            % mean of the folded normal.
            pls_new = pls_new / sqrt(2/pi);
        case 4
            pls_new = abs(s_perm.zsc_Cave);
            if perm2pls_dophi
                pls_new = abs(s_perm.zsc_phi);
            end
            
            if split_plusminus > 0
                sig = get_sig(pvals,perm2pls_do_bh,alpha,bad_any);
                renormalize_on = 0;
                [pls_new] = do_splitplusminus(pls_new,sig,split_plusminus,s_perm,perm2pls_dophi,renormalize_on);
                
            end
    end
    
    if allow_signed
        
        pls_sgn = get_sgn(s_perm,perm2pls_dophi);
        pls_new = pls_new .* pls_sgn;
        
    end
    
end

function pls_new=fix_infs(pls_new)          % Subs some cheesy estimates in place of nans
    sz=size(pls_new);
    temp = pls_new;
    temp(isinf(pls_new)) = 0;
    m = max(temp,[],2);
    m = repmat(m,[1,sz(2),1]);
    m = m*1.2;                          % I estimate 
    pls_new(isinf(pls_new)) = m(isinf(pls_new));

    gaussian_test_inv=0;
    
    if gaussian_test_inv
        % Test for gaussian
        i=10;j=2;
        pls_new_symm = pls_new .* sign(randn(size(pls_new)));       % Simulate making the distribution symmetric about 0 to test for gaussian
        figure; qqplot_with_stats(squeeze(pls_new_symm(i,:,j)));
        figure; hist(squeeze(pls_new_symm(i,:,j)),20);
    end
end


function pls_sgn = get_sgn(s_perm,do_phi)
        if ~do_phi
            dCave = s_perm.Cave1-s_perm.Cave2;
        else
            dCave = s_perm.phi1-s_perm.phi2;
        end
        if do_phi; dCave = get_circ_diff(dCave); end

        pls_sgn = sign(dCave);
end


function sig = get_sig(pvals,perm2pls_do_bh,alpha,bad_any)
    % Calculate significant cells from permutation test
    sz = size(pvals);
    if length(sz) == 2; sz(3)=1; end    % Incase pvals has only 1 ctg, add a 1 to sz(3)
    if perm2pls_do_bh

        sig = false(sz);
        for i = 1:sz(3)
            sig(:,:,i) = bh_step(pvals(:,:,i)',alpha, bad_any)';
        end
    else
        sig = pvals < alpha;
    end
    
    % Remove bad cells
    bad_anyXD = repmat(bad_any(:)',[sz(1), 1, sz(3)]);
    sig(bad_anyXD) = false;
end


function [pls_new,sig_pos,sig_neg] = do_splitplusminus(pls_new0,sig,split_plusminus,s_perm,perm2pls_dophi,renormalize_on,bad_any)
    % Pulls out only significant positive, negative, or both values from
    % pls_new.


    if split_plusminus > 0
        
        pls_sgn = get_sgn(s_perm,perm2pls_dophi);
        %sig_pos = sig .* (pls_sgn >= 0);              % Positive sig cells 
        %sig_neg = sig .* (pls_sgn < 0);               % Negative sig cells
        
        pls_pos = pls_new0 .* (pls_sgn >= 0);              % Only take positive - set negatives to 0
        pls_neg = pls_new0 .* (pls_sgn < 0);               % Only take negative - set positives to 0
        pls_sig = pls_new0 .* sig;                   % Only take significant (this should be equal to sig for return_mode 1)
        

        if renormalize_on
            % This doesn't work as I would like for z-score mode.
            error('This doesnt work as I would like for z-score mode. Should re-evaluate before using.');
            % 
            % The formula for mean is 1/N*sum(Xi). N is the total number of
            % data points. However, since we're only considering a subset
            % (sig_pos/sig_neg/sig), we should instead divide by sum(sig).
            % This renormalization scales all values in pls_new by N /
            % sum(sig) such that the mean operation performed in
            % plot_matrix3D_custstruct is 1/sum(sig) * sum(Xi). Likewise
            % should be done for the variance. (Note - error in mean might
            % still be messed up; should perhaps switch to storing "non
            % sig" entries as Infs.
            Ncells = size(pls_new0,2);
            num_pos = sum(sig_pos,2);
            num_neg = sum(sig_neg,2);
            num_both = sum(sig,2);
            
            num_pos = repmat(num_pos,[1,Ncells,1]);
            num_neg = repmat(num_neg,[1,Ncells,1]);
            num_both = repmat(num_both,[1,Ncells,1]);
            
            pls_pos = pls_pos * sum(~bad_any) ./ num_pos;
            pls_neg = pls_neg * sum(~bad_any) ./ num_neg;
            pls_sig = pls_sig * sum(~bad_any) ./ num_both;
            
        end
        
        switch split_plusminus
            case 1  % Both pos and neg, interspersed
                pls_new = cat(3,zeros(size(pls_new0)),zeros(size(pls_new0)));
                pls_new(:,:,1:2:end) = pls_pos;
                pls_new(:,:,2:2:end) = pls_neg;
            case 2  % Only positive
                pls_new = pls_pos;
            case 3  % Only negative
                pls_new = pls_neg;
            case 4  % Only significant
                pls_new = pls_sig;
        end
        
    
    else
        pls_new = pls_new0;
    end

end