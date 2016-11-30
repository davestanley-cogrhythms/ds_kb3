


function [lfp_wind,tw] = spike_triggered_windows(spikes,lfp,dt,window_size)

    plot_debug = 0;
    plot_on = 0;
    
    Ntot = length(spikes(:));
    t=(0:Ntot-1)*dt;
    Ntrials=size(spikes,2);
    N=size(spikes,1);
    t=reshape(t,N,Ntrials);
    
    
    if plot_on
        figure;
        subplot(211); plot(t,(lfp));
        hold on; plot(t,3*std(lfp(:))*(spikes),'.');
        subplot(212);
        [P f] = pwelch(lfp,window_size,round(window_size/10),window_size,1/dt);
        plot(f,P);
    end
    
    % Remove spikes that are too close to boundaries
    boundary_points = 1:N:(Ntot+1);         % We take Ntot+1 in order to get the final boundary point and to keep spacing consistent
    spikes_crop = crop_edge_spikes(spikes(:),window_size,boundary_points);
    
    % Get an array of indices for each spike
    sp_ind = get_spk_indices(spikes_crop(:),dt);
    
    % Get a matrix of window indices around each spike
    winds = get_window_indices(sp_ind,window_size,Ntot);

    if plot_debug
        figure;
        plot(t,spikes);
        hold on; plot(t(sp_ind),spikes(sp_ind),'kx');
        hold on; plot(t(winds(:)),spikes(winds(:)),'r.');
    end
    

    tw = t(winds);
    sp_wind = spikes(winds);
    lfp_wind = lfp(winds);
    
    if plot_debug
        hold on; plot(t(winds(:,[13])),spikes(winds(:,[13])),'g.');
        hold on; plot(tw(:,1),mean(lfp_wind,2),'k','LineWidth',2);
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



function spikes_crop = crop_edge_spikes(spikes,window_size,boundary_points)

% 
%     % % % Testing mode 1
%     boundary_points = 1:10:31;
%     window_size = 5;
%     spikes = ones(30,1);        % Constant spikes!
%     
%     % % % Testing mode 2
% %     boundary_points = 1:10:31;
% %     window_size = 4;
% %     spikes = ones(30,1);        % Constant spikes!
%     
%     % % % %
%     
    plot_on = 0;

    wh=floor(window_size / 2); % half window
    
    if ~mod(window_size,2)      % Even          
        blackout = (-wh+1):(wh-1);
        % Should be no spikes within d-wh+1 and d+wh-1 inclusive
    else                        % Odd
        blackout = (-wh):(wh-1);
        % Should be no spikes within d-wh and d+wh-1 inclusive
    end
    
    N = length(blackout);
    Npoints = length(boundary_points);

    boundary_all = repmat((boundary_points(:))',N,1) + repmat(blackout(:),1,Npoints);
    boundary_all=boundary_all(:);
    spikes_crop = spikes;
    boundary_all = boundary_all( boundary_all >= 1 & boundary_all <= length(spikes));
    spikes_crop(boundary_all(:)) = 0;
    
    
    if plot_on
        figure; plot(spikes,'.');
        hold on; plot(spikes_crop,'k.');
        plot([boundary_points(:)'; boundary_points(:)'],[ones(1,Npoints);zeros(1,Npoints)],'k','LineWidth',2);
        
    end

end

