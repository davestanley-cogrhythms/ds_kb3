

function sp_ind = get_spk_indices(spikes,dt)

    plot_debug = 0;
    colourarr = 'bgrymcbgrymcbgrymcbgrymcbgrymc';
    
    N = length(spikes);
    
    
    
    
    mxs_bin = max(spikes); % max spikes per bin
    Nspikes = sum(spikes); % total number of spikes
    
    clear sp_ind
    for i = 1:mxs_bin
        sp_ind{i} = (find(spikes(:) >= i))';         % sp_ind of the bin containing at least i spikes. (Most bins contain zero).
    end
    
    if plot_debug
        t=(0:N-1)*dt;
        figure;
        plot(t,spikes);
        for i = 1:mxs_bin
            hold on; plot(t(sp_ind{i}),spikes(sp_ind{i}),[colourarr(i) '.']);
        end
    end
    
    % Merge sp_ind into 1 array. If a single bin contains 5 spikes, for
    % example, then this bin will be referenced 5 times.
    sp_ind = cell2mat(sp_ind);
    sp_ind = sort(sp_ind);
    
    
    
    if plot_debug
        figure;
        plot(t,spikes);
        hold on; plot(t(sp_ind),spikes(sp_ind),'r.');
    end
    
    
    
end






