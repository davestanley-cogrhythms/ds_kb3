

function [C,phi,S12,S1,S2,f,confC,phistd,Cerr] = partcoh(varargin)

    % Variable names
    % R is S12 / (S1*S2)
    % suffix all implies all electrodes

    return_full_matrix = 0; % If set to 1, does the same as partcoh_full
                            % Use for testing only, should set to zero
                            % normally
    
    % Pull out variables from varargin
    fname = varargin{1};
    data = varargin{2};
    chosen_dims = varargin{3};
    if nargin > 3
        args = varargin(4:end);
        params = args{1};
    else args = {};
        params = [];
    end
    [tapers,pad,Fs,fpass,err,trialave]=getparams(params);   % A chronux function
    
    % Calculate chosen dimensions
    Nelects = size(data,3);
    Ntrials = size(data,2);
    chosen_dims_logical = build_logical(chosen_dims,Nelects);
    
    % % Calculate Rall
    % Define data groupings & calc Rxy
    datax = data(:,:,chosen_dims(1));
    datay = data(:,:,chosen_dims(2));

    % Estimate dimensions of output
    N=check_consistency(datax,datay);
    nfft=max(2^(nextpow2(N)+pad),N);
    [f,findx]=getfgrid(Fs,nfft,fpass); 
    Nf = length(f(:));
    
    if ~trialave;
        error('Trialave must be set to 1!! Cannot handle multiple trials');
    end

    Rall = zeros(Nf,Nelects,Nelects);   % R is complex C (i.e. normalized cross spectra)

    for i = 1:Nelects
        for j = 1:Nelects
            [~,~,S12,S1,S2] = fname(data(:,:,i),data(:,:,j),args{:});
            Rall(:,i,j) = S12 ./ sqrt(S1.*S2);
        end
    end
  
    % % Do inverse method
    % Permute Rall so frequencies (f) at the back
    Rall2 = permute(Rall,[2,3,1]);
    d = zeros(size(Rall2));
    for k = 1:Nf
        d(:,:,k,:) = inv(Rall2(:,:,k));   % This inverse method is from Biosignal Processing- Principles and Practices
    end
    
    i=chosen_dims(1);
    j=chosen_dims(2);
    dji = squeeze(d(j,i,:));
    dii = squeeze(d(i,i,:));
    djj = squeeze(d(j,j,:));
    Rxy_z = (-1)^(i+j) .* dji ./ sqrt(djj.*dii);


    % Get C and phi
    C=abs(Rxy_z); 
    phi=angle(Rxy_z);

    % Get S's (don't need to worry about dependency here.)
    [C2,~,~,S1,S2,f] = fname(datax,datay,args{:});
    S12 = Rxy_z.*sqrt(S1.*S2);
    
    % Get other outputs
    if nargout >=9 && err == 2; 
        %warning('Not sure these are correct for partial coh.');
        N=check_consistency(datax,datay);
        nfft=max(2^(nextpow2(N)+pad),N);
        [f,findx]=getfgrid(Fs,nfft,fpass); 
        tapers=dpsschk(tapers,N,Fs); % check tapers
        J1=mtfftc(datax,tapers,nfft,Fs);
        J2=mtfftc(datay,tapers,nfft,Fs);
        J1=J1(findx,:,:); J2=J2(findx,:,:);
        [confC,phistd,Cerr]=coherr(C,J1,J2,err,trialave);
    elseif nargout>=7;
        %warning('Not sure these are correct for partial coh.');
        N=check_consistency(datax,datay);
        nfft=max(2^(nextpow2(N)+pad),N);
        [f,findx]=getfgrid(Fs,nfft,fpass); 
        tapers=dpsschk(tapers,N,Fs); % check tapers
        J1=mtfftc(datax,tapers,nfft,Fs);
        J2=mtfftc(datay,tapers,nfft,Fs);
        J1=J1(findx,:,:); J2=J2(findx,:,:);
        [confC,phistd]=coherr(C,J1,J2,err,trialave);
        Cerr = zeros(2,size(C,1));  % Dummy data for Cerr.
        %warning('Make sure this makes sense');
    end
    
    
    
end

