

function [Call_z,phiall_z,S12all,S1allp,S2allp,f,confC,phistd,Cerr] = partcoh_full(varargin)

    % Variable names
    % R is S12 / (S1*S2)
    % suffix all implies all electrodes
    testing_mode = 0;   % Runs a quick non-matrix version to solve one of the matrix elements, for comparison.
    
    % Pull out variables from varargin
    fname = varargin{1};
    data = varargin{2};
    if nargin > 2
        args = varargin(3:end);
        params = args{1};
    else args = {};
        params = [];
    end
    
    
    [tapers,pad,Fs,fpass,err,trialave]=getparams(params);   % A chronux function
    
    % Calculate chosen dimensions
    Nelects = size(data,3);
    
    % % Calculate Rall
    % Define data groupings & calc Rxy
    data1 = data(:,:,1);
    data2 = data(:,:,2);

    % Estimate dimensions of output
    N=check_consistency(data1,data2);
    nfft=max(2^(nextpow2(N)+pad),N);
    [f,findx]=getfgrid(Fs,nfft,fpass); 
    Nf = length(f(:));
    
    if ~trialave;
        error('Trialave must be set to 1!! Cannot handle multiple trials');
    end

    Rall = zeros(Nf,Nelects,Nelects);   % R is complex C (i.e. normalized cross spectra)
    S1all = zeros(Nf,Nelects);
    S2all = S1all;

    
    for i = 1:Nelects
        for j = 1:Nelects
            [~,~,S12,S1,S2] = fname(data(:,:,i),data(:,:,j),args{:});
            Rall(:,i,j) = S12 ./ sqrt(S1.*S2);
            % Save S's also
            if i == j
                S1all(:,i) = S1;
                S2all(:,j) = S2;
            end
        end
    end
    
    % Permute Rall so frequencies (f) at the back
    Rall2 = permute(Rall,[2,3,1]);
    d = zeros(size(Rall2));
    for k = 1:Nf
        d(:,:,k,:) = inv(Rall2(:,:,k));   % This inverse method is from Biosignal Processing- Principles and Practices
    end
    
    % % Return the full matrix of Rall_z (partial)
    
    % First calculate Dji, Dii, Djj, coeff
    Dji = zeros(size(d));
    Dii = Dji; Djj = Dji;
    diags_d = zeros(Nelects,Nf);
    for k = 1:Nf        % Need to iterate along frequencies
        Dji(:,:,k) = d(:,:,k)';
        diags_d(:,k) = diag(d(:,:,k));
        Dii(:,:,k) = repmat(diags_d(:,k)',Nelects,1);
        Djj(:,:,k) = repmat(diags_d(:,k),1,Nelects);
    end
    [J,I] = meshgrid(1:Nelects,1:Nelects);
    IpJ = I+J;
    coeff = (-1*ones(Nelects,Nelects)).^IpJ;
    coeff = repmat(coeff,[1,1,Nf]);

    % Calculate Rall
    Rall_z2 = coeff .* Dji ./ sqrt(Dii.*Djj);

    % Permute so frequencies back at front
    Rall_z = permute(Rall_z2,[3,1,2]);

    % Calculate S12all
    S1allp = repmat(S1all,[1,1,Nelects]);
    S2allp = repmat(S2all,[1,1,Nelects]);
    S2allp = permute(S2allp,[1,3,2]);
    S12all = Rall_z .* sqrt(S1allp .* S2allp);
    
    Call_z = abs(Rall_z);
    phiall_z = angle(Rall_z);
    
    % Fill these out later.
%     sz = size(Call_z);
%     confC = Inf(sz);
%     phistd = Inf(sz);
%     Cerr = [Inf(sz); Inf(sz)];
    
    confC = [];
    phistd = [];
    Cerr = [];
    
    
    if testing_mode
        % Test against non-matrix mode.
        i=1;
        j=5;
        dji = squeeze(d(j,i,:));
        dii = squeeze(d(i,i,:));
        djj = squeeze(d(j,j,:));
        Rxy_z = (-1)^(i+j) .* dji ./ sqrt(djj.*dii);


        % Get C and phi
        C=abs(Rxy_z); 
        phi=angle(Rxy_z);
        
        % Get S's (don't need to worry about dependency here.)
        data1 = data(:,:,i);
        data2 = data(:,:,j);
        [C2,~,~,S1,S2,f] = fname(data1,data2,args{:});
        S12 = Rxy_z.*sqrt(S1.*S2);
        
        figure
        plot(C);
        hold on; plot(Call_z(:,i,j),'r.');
    end
    
end

