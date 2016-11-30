

function [C,phi,S12,S1,S2,f,confC,phistd,Cerr] = partcoh_recurs(varargin)

    global Rall
    do_fast_precalc_mode = 1;
    
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
    chosen_dims_logical = build_logical(chosen_dims,Nelects);
    
    % Calculate Rall
    %tic
    if ~do_fast_precalc_mode 
        % % Simpler code, but slower
        % Define data groupings & calc Rxy
        datax = data(:,:,chosen_dims(1));
        datay = data(:,:,chosen_dims(2));
        dataz = data(:,:,~chosen_dims_logical);
        Rxy = partcoh_helper(fname,datax,datay,dataz,args);

    else
        % % Alternate faster method to calculate Rxy, using precalc
        % Define data groupings & calc Rxy
        datax = data(:,:,chosen_dims(1));
        datay = data(:,:,chosen_dims(2));
        
        % Estimate dimensions of output
        N=check_consistency(datax,datay);
        nfft=max(2^(nextpow2(N)+pad),N);
        [f,findx]=getfgrid(Fs,nfft,fpass); 
        sz = size(f(:));
        if trialave; sz(2) = 1; else sz(2) = size(data,2); end
        
        Rall = zeros(sz(1),sz(2),Nelects,Nelects);
        
        % Precalculate Rall
        for i = 1:Nelects
            for j = 1:Nelects
                [~,~,S12,S1,S2] = fname(data(:,:,i),data(:,:,j),args{:});
                Rall(:,:,i,j) = S12 ./ sqrt(S1.*S2);
            end
        end

        % Calculate Rxy
%         profile on
        Rxy = partcoh_helper_precalc(chosen_dims(1),chosen_dims(2),find(~chosen_dims_logical));
%         profile off
%         profile viewer
        
    end
    %toc
    
    % Get C and phi
    C=abs(Rxy); 
    phi=angle(Rxy);
    
    % Get S's (don't need to worry about dependency here.)
    [C2,~,~,S1,S2,f] = fname(datax,datay,args{:});
    S12 = Rxy.*sqrt(S1.*S2);
    
    
    
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


function Rxy = partcoh_helper(fname,datax,datay,dataz,args)

    % Recursive formula
    if isempty(dataz)
        [C,phi,S12,S1,S2] = fname(datax,datay,args{:});
        Rxy = S12 ./ sqrt(S1.*S2);
        %Rxy = C .* exp(1i*phi);
    else
        dataz_z0 = dataz(:,:,2:end);
        dataz0 = dataz(:,:,1);
        %partcoh_handle = @(x,y,z) partcoh_helper(fname,x,y,z,args);
        Rxy_z0 = partcoh_helper(fname,datax,datay,dataz_z0,args);
        Rxz0_z = partcoh_helper(fname,datax,dataz0,dataz_z0,args);
        Rz0y_z = partcoh_helper(fname,dataz0,datay,dataz_z0,args);  % Note the order here is important; if do yz0, have to put conj in the equation below
        
%         Rxy = (Rxy_z0 - Rxz0_z .* Ryz0_z) ...
%             ./ (sqrt(1-Rxz0_z.*conj(Rxz0_z)) .* sqrt(1-Ryz0_z.*conj(Ryz0_z)));
        
        Rxy = (Rxy_z0 - Rxz0_z.*Rz0y_z) ./ (sqrt((1-abs(Rxz0_z).^2) .* (1-abs(Rz0y_z).^2)));
        
            
    end

end


function Rxy = partcoh_helper_precalc(dimx,dimy,dimsz)
    global Rall

    % Recursive formula
    if isempty(dimsz)
        Rxy = Rall(:,:,dimx,dimy);  % Query appropriate Rxy
    else
        dimz_z0 = dimsz(2:end);
        dimz0 = dimsz(1);
        
        %partcoh_handle = @(x,y,z) partcoh_helper(fname,x,y,z,args);
        Rxy_z0 = partcoh_helper_precalc(dimx,dimy,dimz_z0);
        Rxz0_z = partcoh_helper_precalc(dimx,dimz0,dimz_z0);
        Rz0y_z = partcoh_helper_precalc(dimz0,dimy,dimz_z0);  % Note the order here is important; if do yz0, have to put conj in the equation below
        
%         Rxy = (Rxy_z0 - Rxz0_z .* Ryz0_z) ...
%             ./ (sqrt(1-Rxz0_z.*conj(Rxz0_z)) .* sqrt(1-Ryz0_z.*conj(Ryz0_z)));
        
        Rxy = (Rxy_z0 - Rxz0_z.*Rz0y_z) ./ (sqrt((1-abs(Rxz0_z).^2) .* (1-abs(Rz0y_z).^2)));
        
            
    end

end