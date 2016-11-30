

function [J0,Jparams] = precalc_fft_and_save(data_all0,os,do_spikes,curr_electrode)

    if ~exist('curr_electrode','var'); curr_electrode = 0;
    end
    
    if curr_electrode > 0
        i0j0_split = 1;
    else
        i0j0_split = 0;
    end

    % Unpack os
    fname = os.fname;
    curr_stage = os.curr_stage;
    params = os.params;
    
    % Initialize stuff
    J0 = [];
    Jparams = [];
    
    
    % Bypass FFT calculation if not needed
    data_type = mode2datatype(os.sfc_mode);
    if strcmp(data_type,'FR'); return; end          % Return if just doing FR
    if do_spikes && ~os.needed_spikes; return; end  % Doing spikes and spikes not needed
    if ~do_spikes && ~os.needed_lfp; return; end    % Doing lfp and lfp not needed
    
    
    % Setup for buildstruct
    s = struct;
    recalc = 0;
    reload = 0;
    
    if nargout == 0; reload = -1; end   % Don't load the FFT data if nargout is zero!
    
    
    
    % Setup path name
    [~,fname_prefix] = fileparts(fname);
    tapers_str = ['tap_' num2str(params.tapers(1)) '_' num2str(params.tapers(2))];
    if os.is_spectrogram
        precalc_path = fullfile(getpath('path_precalcJ_curr'),['stage_' num2str(curr_stage)],['Nw_' num2str(os.Nwind)],['fo_' num2str(os.fract_overlap)],tapers_str);
    else
        precalc_path = fullfile(getpath('path_precalcJ_curr'),['stage_' num2str(curr_stage)],'full',tapers_str);
    end
    
    if do_spikes; fname_out = [fname_prefix '_sp'];
    else fname_out = [fname_prefix '_lfp'];
    end
    
    if i0j0_split
        fname_out = [fname_out '_' num2str(curr_electrode)];
    end
    
    
    % Process FFT's
    tstart = tic;
    fieldname = fname_out;
    myfunc = @() precalc_fft (data_all0,os,do_spikes);
    
    %s = buildstruct_buffered(s, precalc_path, fieldname, myfunc, recalc, reload);
    s.(fieldname) = myfunc();
        % Comparison of processing times shows that recalculating vs
        % reloading FFT coefficients takes about the same amount of time, so
        % we'll go the reloading route. Also, save takes forever so it's
        % best to just recalculate.

    
    
    if nargout > 0
        Js = s.(fieldname);
        J0 = Js.Jout;
        telapsed = toc(tstart);

        Jparams = rmfield(Js,'Jout');
        Jparams.telapsed = telapsed;
    end

end    
    