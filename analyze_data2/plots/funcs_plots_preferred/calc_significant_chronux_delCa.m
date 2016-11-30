

function [sig_cells, pvals] = calc_significant_chronux_delCa(f,Cave,Cerr1,Cerr2,alpha,bad_any,freqval,cols,do_bh_stepup)

    if ~exist('do_bh_stepup','var'); do_bh_stepup = 1; end      % Use BH setup up algorithm to enforce fixed false discovery rate, rather than just simply doing repeated t test
    plot_on = 0;
    plot_debug = 0;
    plot_ctg = 1;        % Note: Ctg 11 may be fake data!
    
    % Move to Fisher space
    Cave = fisher_trans(Cave);
    Cerr1 = fisher_trans(Cerr1);
    Cerr2 = fisher_trans(Cerr2);
    
    if plot_debug       % 
        %% Check for symmetry in fisher space
        figure
        for i = 1:size(Cave,2)
            clf;
            plot([Cave(:,i,plot_ctg) Cerr1(:,i,plot_ctg) Cerr2(:,i,plot_ctg)]);
            pause
        end
    end
    
    
    
    % Look only at a specific frequency
    %index = find(f >= freqval, 1, 'first');
    index = find_closest(f,freqval);
    Cave = squeeze(Cave(index,:,:));
    Cerr1 = squeeze(Cerr1(index,:,:));
    Cerr2 = squeeze(Cerr2(index,:,:));
    
    
    % Get standard deviations
    %alpha = 0.05;    % Alpha of confidence intervals
    Cstd = confid_to_stdev(Cave,Cerr2,alpha);
    if plot_debug
        %% Compare Stdev range
        figure
        plot([Cave(:,plot_ctg) Cerr1(:,plot_ctg) Cerr2(:,plot_ctg)]);
        hold on;
        plot([Cave(:,plot_ctg) Cave(:,plot_ctg)-1.96*Cstd(:,plot_ctg) Cave(:,plot_ctg)+1.96*Cstd(:,plot_ctg)],'k')
    end
    
    z=Cave(:,cols(1)) - Cave(:,cols(2));
    zstd = sqrt(Cstd(:,cols(1)).^2 + Cstd(:,cols(2)));
    
    if ~do_bh_stepup
        [sig_cells,pvals] = norminv_test(z,zstd,alpha); sum(sig_cells)
    else
        [sig_cells,pvals,pcrit] = norminv_test_bh(z,zstd,alpha,bad_any);%sum(sig_cells)
    end
    
    %[sig_cells] = confid_test([Cerr1(:,cols(1)) Cerr2(:,cols(1))],[Cerr1(:,cols(2)) Cerr2(:,cols(2))]);sum(sig_cells)    % Tests for non-overlapping confidence intervals
    
    if plot_on
        figure; plot(pvals)
    end

end

function [sig_cells,pvals] = norminv_test(z,zstd,alpha)
    

    t_statistic = abs(z) ./ zstd;                    % According to Mikio, the t statistic is equal to mu-hat/sig-hat and this follows the t-distribution for 2N degrees of freedom (N=2*tapers*trials)
    

    
    
    Nsigs = norminv(1-alpha/2,0,1);
    %sig_cells2 = t_statistic > Nsigs;                % Furthermore, this t-distribution can be approximated as Gaussian when N is large.

    areas = normcdf(t_statistic,0,1) - normcdf(-t_statistic,0,1);
    pvals = 1 - areas;
    sig_cells = pvals < alpha;


end



function sig_cells = confid_test(x1bounds,x2bounds)     % Tests for non-overlapping confidence intervals

    ind_2_greaterthan_1 = x2bounds(:,1) > x1bounds(:,2);    % Non overlapping when X2 is greater than X1
    ind_1_greaterthan_2 = x1bounds(:,1) > x2bounds(:,2);    % Non overlapping when X1 is greater than X2
    
    sig_cells = ind_2_greaterthan_1 | ind_1_greaterthan_2;  % Either of these conditions is ok!

end
