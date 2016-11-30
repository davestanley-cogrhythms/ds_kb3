

function y = spiketimes_to_timeseries(spikes_indices_jitter,N)

    y = zeros(1,N);
    y(spikes_indices_jitter) = 1;
   
end