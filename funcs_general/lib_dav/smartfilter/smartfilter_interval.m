
function [interval centres] = smartfilter_interval (datatimes, data)

    lf_cutoff = 0; % LF Cutoff at 5 Hz
        lf_val = 20;

    raw_interval = [];
    raw_interval = [60 120 180 240 300 360 420 540 600 660 780 900 1500];    % Default electrical AC filtering


    %Use smart filtering to remove excess mechanical noise
    [f, fft_val] = daveFFT(datatimes, data, 1);
    temp = round(length(f)/2); f = f(1:temp); fft_val = fft_val(1:temp);
    
    notch_bandwidth = 0.5;
    window_hertz = 60;  %Hz
    spike_threshold_multiplier = 20;
    start_freq = 300; %Hz --> Start searching at this frequency
    raw_interval = [raw_interval identify_fft_noise(f, fft_val, start_freq, window_hertz, spike_threshold_multiplier)];
    df = f(2) - f(1);
    

    interval = zeros(length(raw_interval), 2);
    for i = 1:length(raw_interval)
%        interval(i,:) =  [max((raw_interval(i) - 2*df), 0) (raw_interval(i) + 2*df)];   %no negative intervals
       interval(i,:) =  [max((raw_interval(i) - notch_bandwidth), 0) (raw_interval(i) + notch_bandwidth)];   %no negative intervals
    end
    

    centres = raw_interval;
    
    % Add low freq cuttoff
    if lf_cutoff
        interval = [0 lf_val; interval];
        centres = [lf_val/2 centres];
    end
    
end


function [f, X] = daveFFT (t_input, x_input, scale_freq)


    % Scale freq is not used - see daveFFT_scale
    if nargin < 3
        scale_freq = 1;
    end

    tstep = t_input(2)-t_input(1);
    min_t = min(t_input);
    max_t = max(t_input);

    t = min_t:tstep:max_t;
    t = (0:length(t_input)-1)*tstep + min_t;
    N = length(t);
    dt = tstep;
    df = 1/(N*dt);
    %t = (0:N-1)*dt;
    f = df * (0:N-1);

    % x = interp1(t_input,x_input,t);
    x = x_input;
    X = 2*pi*fft (x)/N; % This 2*pi scaling produces the correct amplitudes for CTFT
                        % Dividing by N removes the "delta" functions
    X=X';

end


function noise_bands = identify_fft_noise (f, fft, start_freq, window_hertz, spike_threshold_multiplier)

    start_freq; % Only search for noise at frequencies above this one
    window_hertz; % window size measured in measured in hertz
    spike_threshold_multiplier; % Consider it a spike if it is this many times larger than the average surrounding power

    %Convert this to the indicies of f
    df = f(2) - f(1);
    window = round(window_hertz / df);

    %Change FFT coefficients to power spectrum
    powerSpec = (abs(fft).^2);

    noise_bands = [];
    start_loc = find (f >= start_freq, 1, 'first'); % Identify location of starting frequency
    for i = start_loc:window:(length(f)-window)
        meanpower = mean(abs(powerSpec(i:(i+window-1))));  %Approxixmate power of power spectrum for this particular window
        noise_bands = [noise_bands (find(powerSpec(i:(i+window-1)) > (spike_threshold_multiplier*meanpower)) + (i-1))];   %Find if there are any spikes where the power is many times greater than the background noise

    end

    noise_bands = noise_bands * df;

end

