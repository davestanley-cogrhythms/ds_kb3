
function [C_out,phi_out,S12_out,S1_out,S2_out,f,confC]=postcalc_coherency(trials, func_data)


    % Initialize 
    S12_out = [];
    S1_out = [];
    S2_out = [];
    
    % Func_data
    alldata = func_data.alldata;
    os = func_data.os;
    mypairs = func_data.mypairs;

    % Jparams
    f = alldata.Jparams.f;
    confC = alldata.Jparams.confC;
    N = alldata.Jparams.N;
    nfft = alldata.Jparams.nfft;

    % Parmas
    trialave = 1;

    % OS
    params = os.params;
    baseline_subtract = os.baseline_subtract;
    
    % See if need both J1 and J2
    if os.ue_pairs == 6                      % ue_pairs == 6, we're doing PSDs. Hence, only precalculate J1
        do_J2 = 0;
    else
        do_J2 = 1;
    end

    data_type = mode2datatype(os.sfc_mode);
    switch data_type
        case 'SFC'
            field1 = 'spike';
            field2 = 'lfp';
            do_spikes1 = 1;
            do_spikes2 = 0;
        case 'FFC'
            field1 = 'lfp';
            field2 = 'lfp';
            do_spikes1 = 0;
            do_spikes2 = 0;
        case 'FR'
            error('Incomplete');
        case 'PSD'
            field1 = 'lfp';
            do_spikes1 = 0;
            
    end
        
    if ~os.is_spectrogram
        % Extract
        J1 = alldata.(['J0_' field1])(:,:,trials,mypairs(:,1),:);
        

        % For baseline subtraction
        if baseline_subtract

            data1 = alldata.([field1 '_all0']);
            data1 = mean(data1(:,trials,:,:),2);
            s = precalc_fft (data1,os,do_spikes1);
            clear data1

            J1_base = s.Jout;
            J1_base = repmat(J1_base,[1,1,length(trials),1]);
            J1 = J1 - J1_base(:,:,:,mypairs(:,1),:);
            clear J1_base
        end

        if do_J2
            % Extract
            J2 = alldata.(['J0_' field2])(:,:,trials,mypairs(:,2),:);

            % For baseline subtraction
            if baseline_subtract

                data2 = alldata.([field2 '_all0']);
                data2 = mean(data2(:,trials,:,:),2);
                s = precalc_fft (data2,os,do_spikes2);
                clear data2

                J2_base = s.Jout;
                J2_base = repmat(J2_base,[1,1,length(trials),1]);
                J2 = J2 - J2_base(:,:,:,mypairs(:,2),:);
                clear J2_base

            end
        
            % Complete the coherence calculation
            if do_spikes1 || do_spikes2 
                S12=squeeze(mean(conj(J2).*J1,2));
                S1=squeeze(mean(conj(J2).*J2,2));
                S2=squeeze(mean(conj(J1).*J1,2));
                                        % In Chronux, J1 is lfp and J2 is spikes, which
                                        % is the opposite of what we have here. Hence,
                                        % take the conjugate to make things
                                        % match up.
            else
                S12=squeeze(mean(conj(J1).*J2,2));
                S1=squeeze(mean(conj(J1).*J1,2));
                S2=squeeze(mean(conj(J2).*J2,2));
            end
            clear J1 J2

            if trialave; S12=squeeze(mean(S12,2)); S1=squeeze(mean(S1,2)); S2=squeeze(mean(S2,2)); end;

            if strcmp(data_type,'SFC') && os.thinning_mode == 1
                rate0 = alldata.([field1 '_all0']);
                rate1 = squeeze(mean(mean(rate0(:,trials,:),1),2)/get_dt);
                rate2a_old = repmat(rate1',[size(S12,1),1]);
                rate2 = rate1(mypairs(:,1))';
                rate3 = repmat(rate2,[size(S12,1),1]);
%                 rate = repmat()',[size(S12,1),1]);
                %====Aoi code=============
                StarRate = get_rateStar();
                c0 = StarRate./rate3;
                C12=S12./sqrt(S1.*S2)./sqrt((1-c0).*rate3./S2./c0+1);
            else
                C12=S12./sqrt(S1.*S2);
            end

            C=abs(C12); 
            phi=angle(C12);

            C_out = C;
            phi_out = phi;
    %         S12_out = S12;
    %         S1_out = S1;
    %         S2_out = S2;
        else
            % Just doing PSD in this case
            S1=squeeze(mean(conj(J1).*J1,2));
            if trialave; S1=squeeze(mean(S1,2)); end
            C_out = S1;
            phi_out = [];
        end
        
    else
        
        if baseline_subtract
            % Calculate baselines
            data1 = alldata.([field1 '_all0']);
            data1 = mean(data1(:,trials,:,:),2);
            s = precalc_fft (data1,os,do_spikes1);
            clear data1
            J1_base0 = s.Jout;
            clear s
            
            
            if do_J2
                % Calculate baselines
                data2 = alldata.([field2 '_all0']);
                data2 = mean(data2(:,trials,:,:),2);
                s = precalc_fft (data2,os,do_spikes2);
                clear data2
                J2_base0 = s.Jout;
                clear s
            end
        end
        
        chunk_size = 3;
        Npairs = size(mypairs,1);
        Nchunks = ceil(Npairs / chunk_size);
        
        for ch = 1:Nchunks
            j = chunk_size*(ch-1)+1:min(chunk_size*ch,Npairs);
            
            % Extract
            J1 = alldata.(['J0_' field1])(:,:,trials,:,mypairs(j,1));

            % For baseline subtraction
            if baseline_subtract
                J1_base = repmat(J1_base0(:,:,:,:,mypairs(j,1)),[1,1,length(trials)]);
                J1 = J1 - J1_base;
                clear J1_base_temp
            end
            
            if do_J2
                % Extract
                J2 = alldata.(['J0_' field2])(:,:,trials,:,mypairs(j,2));

                % For baseline subtraction
                if baseline_subtract
                    J2_base = repmat(J2_base0(:,:,:,:,mypairs(j,2)),[1,1,length(trials)]);
                    J2 = J2 - J2_base;
                    clear J2_base_temp
                end

                % Complete the coherence calculation
                if do_spikes1 || do_spikes2 
                    S12=squeeze(mean(conj(J2).*J1,2));
                    S1=squeeze(mean(conj(J2).*J2,2));
                    S2=squeeze(mean(conj(J1).*J1,2));
                                            % In Chronux, J1 is lfp and J2 is spikes, which
                                            % is the opposite of what we have here. Hence,
                                            % take the conjugate to make things
                                            % match up.
                else
                    S12=squeeze(mean(conj(J1).*J2,2));
                    S1=squeeze(mean(conj(J1).*J1,2));
                    S2=squeeze(mean(conj(J2).*J2,2));
                end
                clear J1 J2

                if trialave; S12=squeeze(mean(S12,2)); S1=squeeze(mean(S1,2)); S2=squeeze(mean(S2,2)); end;

                if strcmp(data_type,'SFC') && os.thinning_mode == 1
                    rate = alldata.([field1 '_all0']);
                    rate = squeeze(mean(mean(rate(:,trials,:,mypairs(j,1)),1),2)/get_dt);
                    rate = repmat(rate,[1,1,size(S12,1)]); rate = permute(rate,[3,1,2]);
                    %====Aoi code=============
                    StarRate = get_rateStar();
                    c0 = StarRate./rate;
                    C12=S12./sqrt(S1.*S2)./sqrt((1-c0).*rate./S2./c0+1);
                else
                    C12=S12./sqrt(S1.*S2);
                end

                C=abs(C12);
                phi=angle(C12);
            else
                % Just doing PSD in this case
                S1=squeeze(mean(conj(J1).*J1,2));
                if trialave; S1=squeeze(mean(S1,2)); end
                C = S1;
            end
            
            
            
%             C_curr = C_curr;
%             phi_curr = phi_curr;
            
            % Preallocate zeros
            if ~exist('C_out','var'); C_out = zeros([size(C,1),size(C,2),size(mypairs,1)]); end
            if do_J2
                if ~exist('phi_out','var');
                    phi_out = C_out;    % This used to be C; I changed it to C_out but didn't test; might cause error!
                end
            end
%                 S12_out = C;
%                 S1_out = C;
%                 S2_out = C;

            
            C_out(:,:,j) = C;
            if do_J2; phi_out(:,:,j) = phi; end
%             S12_out(:,j,:) = S12;
%             S1_out(:,j,:) = S1;
%             S2_out(:,j,:) = S2;

            
        end
        
%         C_out = permute(C_out,[1,3,2]);
%         phi_out = permute(phi_out,[1,3,2]);
        
        sz_new = [size(C_out,1)*size(C_out,2),size(C_out,3)];
        C_out = reshape(C_out,sz_new);
        if do_J2; phi_out = reshape(phi_out,sz_new); 
        else phi_out = [];
        end
        
    end

end