
%% Generate 8-80Hz coupled signals
    
    Ntrials = 10;
    trial_dur = 1;
    dt=1e-3;
    Fs = 1/dt;
    N = Ntrials*trial_dur;
%     N=10;
    
    freq1=2.5*pi;
    freq2=80;
%     freq1=10;
%     freq2=18;
    mysig = gen_simsig(freq1,freq2,Fs,N,1,10);
    t = 1:length(mysig); t=t/Fs;
    mysig = mysig + 80*sin(2*pi*freq1*t)';
%     mysig = reshape(mysig,trial_dur*Fs,Ntrials);
    


    % Assign input data
    sig_pac = mysig;
    %sig_pac = currlfp1_ctg(:); sig_pac = sig_pac(1:10000);
    %sig_pac = currlfp1_ctg(:,1:100);
    
%     figure; plot(sig_pac);
    
% Run analysis

    sig_pac = sig_pac;
    Fs = 1/dt;
    measure = 'mi';

    % Default values
    sig_mod = sig_pac;
    ph_freq_vec = 2:3:30;
    amp_freq_vec = 10:10:140;
%     ph_freq_vec = [5 8];       % In good region (should be sig)
%     amp_freq_vec = [70 80];
%     ph_freq_vec = [5 8];       % In bad region (should be not significant)
%     amp_freq_vec = [30 40];
%   
%     ph_freq_vec = 2:4:14;
%     ph_freq_vec = [6 10];
%     amp_freq_vec = [75 85];

    plt = 'y';
    waitbar = 1;
    width = 7;
    nfft = ceil(Fs/(diff(ph_freq_vec(1:2))));
    num_shf = 0;
    alpha = 0.05;
    dataname = '';
    sig_pac_name = '';
    sig_mod_name = '';
    
    
    [pacmat, freqvec_ph, freqvec_amp,pmat] = find_pac_shf(sig_pac, Fs, measure, sig_mod, ph_freq_vec, amp_freq_vec, plt, waitbar, width, nfft, num_shf, alpha, dataname, sig_pac_name, sig_mod_name);
%     caxis([0 1])
    



%% Generate alpha-beta coupling
    
    Ntrials = 10;
    trial_dur = 0.8;
    dt=1e-3;
    Fs = 1/dt;
    N = Ntrials*trial_dur;
%     N=10;
    
    freq1=pi*3;
    freq2=80;
%     freq1=10;
%     freq2=18;
    mysig = gen_simsig(freq1,freq2,Fs,N,1,10);
    t = 1:length(mysig); t=t/Fs;
    mysig = mysig + 80*sin(2*pi*freq1*t)';
    mysig = reshape(mysig,trial_dur*Fs,Ntrials);
    


    % Assign input data
    sig_pac = mysig;
    %sig_pac = currlfp1_ctg(:); sig_pac = sig_pac(1:10000);
    %sig_pac = currlfp1_ctg(:,1:100);
    
    
    if size(sig_pac,2) > 1
        figure; plot(sig_pac(:,1:1));
    else
        figure; plot(sig_pac);
    end
    xlabel('time (ms)'); ylabel('mV');
    

% Run analysis

    sig_pac = sig_pac;
    Fs = 1/dt;
    measure = 'mi';

    % Default values
    sig_mod = sig_pac;
    ph_freq_vec = 5:1:25;
    amp_freq_vec = 5:1:90;
%     ph_freq_vec = [5 8];       % In good region (should be sig)
%     amp_freq_vec = [70 80];
%     ph_freq_vec = [5 8];       % In bad region (should be not significant)
%     amp_freq_vec = [30 40];
%   
%     ph_freq_vec = 2:4:14;
%     ph_freq_vec = [6 10];
%     amp_freq_vec = [75 85];

    plt = 'y';
    waitbar = 1;
    width = 7;
    nfft = ceil(Fs/(diff(ph_freq_vec(1:2))));
    num_shf = 0;
    alpha = 0.05;
    dataname = '';
    sig_pac_name = '';
    sig_mod_name = '';
    
    
    [pacmat, freqvec_ph, freqvec_amp] = find_pac_shf(sig_pac, Fs, measure, sig_mod, ph_freq_vec, amp_freq_vec, plt, waitbar, width, nfft, num_shf, alpha, dataname, sig_pac_name, sig_mod_name);
    ind  = find(freqvec_amp > 20, 1, 'first');
    cmax = max(max(pacmat(ind:end,:)))
    caxis([0 300])
    
    
%% Original defaults
% Necessary inputs values
    sig_pac = sig_pac(:,:);
    Fs = 1/dt;
    measure = 'mi';

% Default values
    sig_mod = sig_pac;
    ph_freq_vec = 1:5:101;
    amp_freq_vec = 1:5:101;
    plt = 'y';
    waitbar = 1;
    width = 7;
    nfft = ceil(Fs/(diff(ph_freq_vec(1:2))));
    num_shf = 0;
    alpha = 0.05;
    dataname = '';
    sig_pac_name = '';
    sig_mod_name = '';
    
    
    [pacmat, freqvec_ph, freqvec_amp] = find_pac_shf (sig_pac, Fs, measure, ...
    sig_mod, ph_freq_vec, amp_freq_vec, plt, waitbar, width, nfft, num_shf, alpha,...
    dataname, sig_pac_name, sig_mod_name);
