


function [sfc f lfp_wind] = sfc_fries2D(spikes,lfp,dt,window_size)

    % Window size must be in indices.

    
    colourarr = 'bgrymcbgrymcbgrymcbgrymcbgrymc';
    

    [lfp_wind,tw] = spike_triggered_windows(spikes,lfp,dt,window_size);

    % calculate spike field coherence based on all the windows
    if ~isempty(lfp_wind)
        [sfc f] = calc_sfc_core(lfp_wind,dt);
    else
        sfc = zeros(floor(size(lfp_wind,1)/2)+1,1)*NaN;
        f = sfc;
    end
    

    
    
end


function [sfc f] = calc_sfc_core(lfp_wind,dt)
    plot_debug = 0;
    fs = 1/dt;
    
    lfp_wind = zscore(lfp_wind);
    Nwind = size(lfp_wind,2);
    
    
    lfp_mu = mean(lfp_wind,2);
    lfp_std = std(lfp_wind,[],2);
    
    % % Calculate numerator
    [Pnum] = psd_wrapper_fries(lfp_mu,fs);
    
    % % Calculate denominator
%     tic
    [Pmat f] = psd_wrapper_fries(lfp_wind,fs);
%     toc
    
    psd_mu = mean(Pmat,2);
    psd_std = std(Pmat,[],2);
    
    sfc = Pnum ./ psd_mu;
    
    if plot_debug
        show_raw=1;
        max_number_to_show = 100;
            if max_number_to_show ==0; show_raw=0; end
            
        N = size(lfp_wind,1);
        t=(0:N-1)*dt;
        N_showing= min([size(lfp_wind,2) max_number_to_show]);
        figure('Position',[ 28         355        1028         576]);
        subplot(221);
        if show_raw; h1 = plot(t,lfp_wind(:,1:N_showing)); end
        %hold on; errorbar(t(:),lfp_mu,lfp_std,'k','LineWidth',2);
        hold on; h2 = plot(t(:),lfp_mu,'k','LineWidth',2);
        if show_raw; h = [h1(1) h2(1)]; legend(h,'Single','Mean'); end
        title('Event triggered');

        subplot(222);
        hold on; plot(f,Pnum,'b'); xlim([0 100]);
        %hold on; errorbar(f(:),psd_mu,psd_std,'r','LineWidth',2);
%         hold on; plot(f(:),psd_mu,'r','LineWidth',2);
        legend('Numerator');
        %title('Numerator vs denominator');
        subplot(223);
        if show_raw; hold on; h1 = plot(f,Pmat(:,1:N_showing),'b'); end
        %hold on; errorbar(f(:),psd_mu,psd_std,'r','LineWidth',2);
        hold on; h2 = plot(f(:),psd_mu,'r','LineWidth',2); xlim([0 100]);
        if show_raw; h = [h1(1) h2(1)]; legend(h,'Mean denominator','All denomin.'); end
        subplot(224);
        hold on; plot(f,sfc); xlim([0 100]);
        legend('Spike field coherence');
    end



end




function [P, f] = psd_wrapper_fries(lfp,fs, mode)

    %if ~exist('N','var'); N = length(lfp); end
    if ~exist('fs','var'); fs = 1; end
    if ~exist('mode','var'); mode = 2; end          % Use Chronux mode by defualt
    
    
    if isvector(lfp)
        lfp=lfp(:);
    end
    
    N = size(lfp,1);
    Nwinds = size(lfp,2);
    
    switch mode
        case 1
            [Ptemp f] = pwelch(lfp(:,1),N,1,N,fs);
            ln = length(Ptemp); clear Ptemp     % Test run to get length
            P = zeros(ln,Nwinds);
            % Pwelch with only 1 window. We'll be averaging over trials anyways.
            for j = 1:Nwinds
                [P(:,j) f] = pwelch(lfp(:,j),N,1,N,fs); 
            end
            %P = mean(P,2);
    
        case 2
            params.Fs = fs;
            %params.trialave = 0;
            params.tapers = get_tapers_fries;
            %tic
            [P f] = mtspectrumc(lfp,params); 
            %toc
    end

end




function tapers = get_tapers_fries
    % Returns the default value of tapers that are used throughout
    % the code.
    tapers = [3 5];
    %tapers = [5 9];
end

