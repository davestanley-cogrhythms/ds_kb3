
function s = precalc_fft (data,os,do_spikes)
    %% precalc_fft
    Jout=[];
    
    % Determine whether running on GPU
    S=whos('data');
    is_gpu = strcmp(S.class,'gpuArray');
    

    % Ue pairs
    params = os.params;
    
    % Nelects, Ntimebins
    Nunits = size(data,3);
    Ntimebins = size(data,4);
    
    
    if ~isempty(data)
        [sz] = mtfft_outsize(data(:,:,1),params);
        
        if is_gpu; Jout = zeros([sz, Nunits, Ntimebins],'gpuArray');
        else Jout = zeros([sz, Nunits, Ntimebins]);
        end
        
        % Preconvert tapers
        [tapers,pad,Fs,fpass,err,trialave]=getparams(params);
        N=size(data,1);
        params.tapers=dpsschk(tapers,N,Fs); % check tapers
        if is_gpu params.tapers = gpuArray(params.tapers);
        end
        
        for i = 1:Nunits
            for j = 1:Ntimebins
                if do_spikes
                    [J_temp,f,confC,N,nfft] = mtfftpb_wrapper(data(:,:,i,j),params);
                else
                    [J_temp,f,confC,N,nfft] = mtfftc_wrapper(data(:,:,i,j),params);
                end
                Jout(:,:,:,i,j) = J_temp;
            end
        end
        
    end
    
    s.Jout = Jout;
    s.f = f;
    s.confC = confC;
    s.N = N;
    s.nfft = nfft;
    
    clear J_temp
    
end



function [J1,f,confC,N,nfft] = mtfftpb_wrapper(data1,params)
    %% mtfftpb_wrapper
    % Wrapper function to implement just mt fourier transform, with
    % necessary precalculations.

    data1=change_row_to_column(data1);
    [tapers,pad,Fs,fpass,err,trialave]=getparams(params);
    
    N=size(data1,1);
    nfft=max(2^(nextpow2(N)+pad),N);
    [f,findx]=getfgrid(Fs,nfft,fpass); 
    tapers=dpsschk(tapers,N,Fs); % check tapers
    J1=mtfftpb(data1,tapers,nfft);
    J1=J1(findx,:,:);
    
    p=err(2);
    pp=1-p/2;
    
    [nf,K,Ch]=size(J1);
    dim=K*Ch;
    dof=2*dim;
    
    if dof <= 2
       confC = 1;
    else     
       df = 1./((dof/2)-1);
       confC = sqrt(1 - p.^df);
    end;

end


function [J1,f,confC,N,nfft] = mtfftc_wrapper(data1,params)
    %% mtfftc_wrapper
    % Wrapper function to implement just mt fourier transform, with
    % necessary precalculations.

    data1=change_row_to_column(data1);
    [tapers,pad,Fs,fpass,err,trialave]=getparams(params);
    
    N=size(data1,1);
    nfft=max(2^(nextpow2(N)+pad),N);
    [f,findx]=getfgrid(Fs,nfft,fpass); 
    tapers=dpsschk(tapers,N,Fs); % check tapers
    J1=mtfftc(data1,tapers,nfft,Fs);
    J1=J1(findx,:,:);
    
    p=err(2);
    pp=1-p/2;
    
    [nf,K,Ch]=size(J1);
    dim=K*Ch;
    dof=2*dim;
    
    if dof <= 2
       confC = 1;
    else     
       df = 1./((dof/2)-1);
       confC = sqrt(1 - p.^df);
    end;

end


function [sz] = mtfft_outsize(data1,params)
    %% mtfft_outsize
    % Wrapper function to implement just mt fourier transform, with
    % necessary precalculations.

    data1=change_row_to_column(data1);
    [tapers,pad,Fs,fpass,err,trialave]=getparams(params);
    
    N=size(data1,1);
    nfft=max(2^(nextpow2(N)+pad),N);
    [f,findx]=getfgrid(Fs,nfft,fpass);
    sz = [length(f),params.tapers(2),size(data1,2)];
end