

function [sig_cells, z, zstd, t_statistic] = calc_significant_chronux_phi_bh_stepup(f,phiave,phistd,alpha,bad_any,freqval,cols,do_bh_stepup)

    if ~exist('do_bh_stepup','var'); do_bh_stepup = 1; end      % Use BH setup up algorithm to enforce fixed false discovery rate, rather than just simply doing repeated t test
    plot_on = 0;

    phiave(phiave<0) = phiave(phiave<0) + 2*pi;     % Shift from being -pi to pi -> 0 to 2pi, since things tend to cluster around 180
   
    %index = find(f >= freqval, 1, 'first');
    index = find_closest(f,freqval);
    phiave2 = squeeze(phiave(index,:,:));
    phistd2 = squeeze(phistd(index,:,:));
    
    if plot_on
        figure; imagesc(phistd(:,:,17),[0 20]); title('Phi std, column 17'); xlabel('Cell #');ylabel('Frequency');
        figure; bar(mean(phistd2)); title('Mean phistd across cels')


        ztemp = phiave(:,:,cols(1)) - phiave(:,:,cols(2));
        ind=ztemp>pi; ztemp(ind)=ztemp(ind)-2*pi; ind=ztemp<-pi; ztemp(ind)=ztemp(ind)+2*pi;        % Bounds of maximum shift should be -pi to pi
        ztemp = abs(ztemp);
        figure; imagesc(ztemp); title('Phase shift (z), across selected columns'); xlabel('Cell #'); ylabel('Frequency');
    end
    
    
    z = phiave2(:,cols(1)) - phiave2(:,cols(2));
    ind=z>pi; z(ind)=z(ind)-2*pi; ind=z<-pi; z(ind)=z(ind)+2*pi;        % Bounds of maximum shift should be -pi to pi
    z = abs(z);

    
    zstd = sqrt(phistd2(:,cols(1)).^2 + phistd2(:,cols(2)).^2);
    
    
    
    %sig_cells = abs(z) > (Nsigs * zstd + 0);        % MY way
    
    t_statistic = abs(z) ./ zstd;                    % According to Mikio, the t statistic is equal to mu-hat/sig-hat and this follows the t-distribution for 2N degrees of freedom (N=2*tapers*trials)
    
    if do_bh_stepup
        areas = normcdf(t_statistic,0,1) - normcdf(-t_statistic,0,1);
        pvals = 1 - areas;
        
        %[sig_cells1 p01] = bh_step1D(pvals,alpha);    % My code
        %sig_cells2 = fdr_bh(pvals,alpha);  % Mathworks code - they seem to give the same result
        
        %[~,p0] = bh_step1D(pvals(~bad_any),alpha);    % Account for bad cells
        %sig_cells3 = pvals <= p0;
        
        [sig_cells,p0] = bh_step1D(pvals,alpha,bad_any);    % Account for bad cells
        
    else
        Nsigs = norminv(1-alpha/2,0,1);
        sig_cells = t_statistic > Nsigs;                % Furthermore, this t-distribution can be approximated as Gaussian when N is large.
    end
    
end

