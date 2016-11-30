
%% Main function
function [spike_all0,lfp_all0,T,J0_spike,J0_lfp,Jparams] = preload_and_save_tiled(os,md)
    % If doing spectrogram, don't have enough memory to store everything in
    % RAM. Therefore, save to file and load during execution.
    
    % Unpack os
    fname = os.fname;
    curr_stage = os.curr_stage;
    Nunits = length(md.unit_names);
    Nelects = length(md.lfp_names);
    
    [~,fname_prefix] = fileparts(fname);

    % Initialize outputs
    spike_all0 = [];
    lfp_all0 = [];
    J0_spike = [];
    J0_lfp = [];
    Jparams = [];
    
    % If we are in spectrogram mode, only calculate the data if cache files
    % don't exist yet
    s = struct;
    recalc = 0;
    reload = 0;

    precalc_path = fullfile(getpath('path_precalcD_curr'),['stage_' num2str(curr_stage)],['Nw_' num2str(os.Nwind)],['fo_' num2str(os.fract_overlap)]);
    %precalc_path = strrep(precalc_path,'.','_');
    fname_spike = [fname_prefix '_spk'];
    fname_lfp = [fname_prefix '_lfp'];

    
    % Spike files
    if os.needed_spikes
        % Load all spikes
        data = get_spike_ts_all(Nunits,curr_stage,md);

        for i = 1:Nunits
            % Do timeseries
            fieldname = [fname_spike '_' num2str(i)];
            myfunc = @() load_spikes_tiled(data(:,:,i),os);
            s.(fieldname)=myfunc();
            %s = buildstruct_buffered(s, precalc_path, fieldname, myfunc, recalc, reload);
            
            % Add to output array
            out_temp = s.(fieldname).spike_all0;
            if isempty(spike_all0); spike_all0 = zeros([size(out_temp,1),size(out_temp,2),size(out_temp,3),Nunits]); end
            spike_all0(:,:,:,i) = out_temp;
            clear out_temp
            
            % Do FFT coefficients
            curr_unit = i;
            do_spikes = 1; [out_temp, Jparams] = precalc_fft_and_save(s.(fieldname).spike_all0,os,do_spikes,curr_unit);
            
            % Add to output array
            if isempty(J0_spike); J0_spike = zeros([size(out_temp,1) size(out_temp,2), size(out_temp,3), size(out_temp,4),Nunits]); end
            J0_spike(:,:,:,:,i) = out_temp;
            clear out_temp
        end
        T = s.(fieldname).T;
    end
    

    % LFP files
    if os.needed_lfp
        % Load all lfp
        data = load_lfp(fname,curr_stage,md);

        for i = 1:Nelects
            % Do timeseries
            fieldname = [fname_lfp '_' num2str(i)];
            myfunc = @() load_lfp_tiled(data(:,:,i),os);

            %s = buildstruct_buffered(s, precalc_path, fieldname, myfunc, recalc, reload);   
            s.(fieldname)=myfunc();
                % Calculating directly is now faster because (1) I've
                % optimized this code and (2) it avoids loading and saving
                % to the disk.

            % Add to output array
            out_temp = s.(fieldname).lfp_all0;
            if isempty(lfp_all0); lfp_all0 = zeros([size(out_temp,1),size(out_temp,2),size(out_temp,3),Nelects]); end
            lfp_all0(:,:,:,i) = out_temp;
            clear out_temp
            
            % Do FFT coefficients
            curr_electrode = i;
            do_spikes = 0; [out_temp, Jparams] = precalc_fft_and_save(s.(fieldname).lfp_all0,os,do_spikes,curr_electrode);
            
            % Add to output array
            if isempty(J0_lfp); J0_lfp = zeros([size(out_temp,1) size(out_temp,2), size(out_temp,3), size(out_temp,4),Nelects]); end
            J0_lfp(:,:,:,:,i) = out_temp;
            clear out_temp
            
            % Add to output array
        end
        T = s.(fieldname).T;
    end

end

%% Tiling functions
function [s] = load_spikes_tiled(data,os)
    
    %data = data(1:800,1:92,:); warning('shortening data');
    [data,T] = tile_data_spectrogram(data,os);
    s.spike_all0 = data;
    s.T = T;
end


function [s] = load_lfp_tiled(data,os)
    
    %data = data(1:300,1:92,:); warning('shortening data');
    [data,T] = tile_data_spectrogram(data,os);
    s.lfp_all0 = data;
    s.T = T;
end


function [data_tiled,T] = tile_data_spectrogram(data,os)
    % function [spike_all0,lfp_all0] = tile_data_spectrogram(spike_all0,lfp_all0,os);
    % Supporting function for implementing tiling
    
    plot_check = 0;
    
    ide = @(X) X; % Identity function
    

    Nelects = size(data,3);

    % Do first one
    i=1;
    dat_curr = data(:,:,1);
    [out , T] = func2spectrogram(dat_curr,'fname',ide,'Nwind',os.Nwind,'fract_overlap',os.fract_overlap,'trialave',0);
    sz = size(out);

    data_tiled = zeros_gpuwrapper([sz,Nelects],data);
    data_tiled(:,:,:,1) = out;

    % Do the rest
    for i = 2:Nelects
        dat_curr = data(:,:,i);
        [data_tiled(:,:,:,i) , T] = func2spectrogram(dat_curr,'fname',ide,'Nwind',os.Nwind,'fract_overlap',os.fract_overlap,'trialave',0);
    end

    %data_tiled = permute(data_tiled,[1,2,4,3]);        % For purposes of
                                                        % squeeze screwing up when we have singleton units or electrodes, we
                                                        % want electrodes at the end

    if plot_check
        figure; plot(data(:,1,1));
        temp=squeeze(data_tiled(:,1,1,:));
        hold on; plot(temp(:),'r.');
        figure; imagesc(temp);
    end
    clear Nelects

    %data_tiled = single(data_tiled);
    
end




