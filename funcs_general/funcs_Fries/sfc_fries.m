


function [sfc f] = sfc_fries(v,dt,window_size)

    % Window size must be in indices.

    plot_debug = 0;
    plot_on = 0;
    colourarr = 'bgrymcbgrymcbgrymcbgrymcbgrymc';
    
    spikes = v(:,1);
    lfp = v(:,2);
    N = length(spikes);
    t=(0:N-1)*dt;
    
    if plot_on
        figure;
        subplot(211); plot(t,zscore(lfp));
        hold on; plot(t,zscore(spikes),'r');
        subplot(212);
        [P f] = pwelch(lfp,window_size,round(window_size/10),window_size,1/dt);
        plot(f,P);
    end
    
    % Get an array of indices for each spike
    sp_ind = get_spk_indices(spikes,dt);
    
    % Get a matrix of window indices around each spike
    winds = get_window_indices(sp_ind,window_size,N);

    if plot_debug
        figure;
        plot(t,spikes);
        hold on; plot(t(winds(:)),spikes(winds(:)),'r.');
    end
    

    tw = t(winds);
    sp_wind = spikes(winds);
    lfp_wind = lfp(winds);
    
    if plot_debug
        hold on; plot(t(winds(:,[13])),spikes(winds(:,[13])),'g.');
    end
    


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
    
    Nwind = size(lfp_wind,2);
    
    
    lfp_mu = mean(lfp_wind,2);
    lfp_std = std(lfp_wind,[],2);
    
    % % Calculate numerator
    [Pnum] = psd_wrapper_fries(lfp_mu,fs);
    
    % % Calculate denominator
    [Ptemp] = psd_wrapper_fries(lfp_wind(:,1),fs);
    ln = length(Ptemp); clear Ptemp
    
    Pmat = zeros(ln,Nwind);
    for i = 1:Nwind
        [Pmat(:,i) f] = psd_wrapper_fries(lfp_wind(:,i),fs);
    end
    
    psd_mu = mean(Pmat,2);
    psd_std = std(Pmat,[],2);
    
    sfc = Pnum ./ psd_mu;
    
    if plot_debug
        N = size(lfp_wind,1);
        t=(0:N-1)*dt;
        N_showing= min([size(lfp_wind,2) 100]);
        figure;
        subplot(221);
        h1 = plot(t,lfp_wind(:,1:N_showing));
        %hold on; errorbar(t(:),lfp_mu,lfp_std,'k','LineWidth',2);
        hold on; h2 = plot(t(:),lfp_mu,'k','LineWidth',2);
        h = [h1(1) h2(1)]; legend(h,'Single','Mean')
        title('Event triggered');

        subplot(222);
        hold on; plot(f,Pnum,'b');
        %hold on; errorbar(f(:),psd_mu,psd_std,'r','LineWidth',2);
        hold on; plot(f(:),psd_mu,'r','LineWidth',2);
        legend('Numerator','Denominator');
        title('Numerator vs denominator');
        subplot(223);
        hold on; h1 = plot(f,Pmat(:,1:N_showing),'b');
        %hold on; errorbar(f(:),psd_mu,psd_std,'r','LineWidth',2);
        hold on; h2 = plot(f(:),psd_mu,'r','LineWidth',2);
        h = [h1(1) h2(1)]; legend(h,'Mean denominator','All denomin.');
        subplot(224);
        hold on; plot(f,sfc);
        legend('Spike field coherence');
    end



end










function winds = get_window_indices(sp_ind,window_size,N)


    %window_size = 10;
    wh=floor(window_size / 2); % half window
    
    % Make sure when taking windows we don't exceed bounds of signal
    if ~mod(window_size,2)     % Even          
        sp_ind = sp_ind(sp_ind - wh  > 0);
        sp_ind = sp_ind(sp_ind + wh - 1 <= N);
    else                    % Odd
        sp_ind = sp_ind(sp_ind - wh > 0);
        sp_ind = sp_ind(sp_ind + wh <= N);
    end
    
    Ngs = length(sp_ind);
    winds = (1:window_size) - wh - 1;
    winds = repmat(winds(:),1,Ngs);
    winds = winds + repmat((sp_ind(:))',window_size,1);
    
    

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
            params.tapers = get_tapers;
            %tic
            [P f] = mtspectrumc(lfp,params); 
            %toc
    end

end

