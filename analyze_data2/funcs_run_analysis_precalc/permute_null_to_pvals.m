
function sfc = permute_null_to_pvals(sfc)
    
    for i = 1:length(sfc)

        s = sfc{i};
        if isfield(s,'Cave1');
            s.pvals_Cave = calcp(s.Cave1,s.Cave2,s.Cavep1,s.Cavep2,0);
            if isfield(s,'phi1');
                s.pvals_phi = calcp(s.phi1,s.phi2,s.phip1,s.phip2,1);
            end

            s.zsc_Cave = calc_normalized_dCave(s.Cave1,s.Cave2,s.Cavep1,s.Cavep2,0);
            if isfield(s,'phi1');
                s.zsc_phi = calc_normalized_dCave(s.phi1,s.phi2,s.phip1,s.phip2,1);
            end


            if isfield(s,'phi1');
                s = rmfield(s,{'Cavep1','Cavep2','phip1','phip2'});
            else
                s = rmfield(s,{'Cavep1','Cavep2'});
            end
        end

        sfc{i} = s;
    end
    
end


function [pvals, critval]= calcp(Cave1,Cave2,Cavep1,Cavep2, do_phi)

    % Get diff
    dCave = Cave1 - Cave2;
    if do_phi; dCave = get_min_phi_diff(dCave); end
    dCave = abs(dCave);

    % Git diff null
    dCave_null = Cavep1 - Cavep2;
    if do_phi; dCave_null = get_min_phi_diff(dCave_null); end
    sz = size(dCave_null);
    
    % Centre both values about zero (new code)
    null_means = mean(dCave_null,3);
    dCave = dCave - null_means;
    dCave_null = dCave_null - repmat(null_means,[1,1,sz(3)]);

    
    %pvals = sum(dCave_null > abs(repmat(dCave,[1,1,sz(3)])),3) ./ (sz(3)*ones([sz(1:2)])); % One tailed
    pvals = sum( (dCave_null > abs(repmat(dCave,[1,1,sz(3)]))) | (dCave_null < -abs(repmat(dCave,[1,1,sz(3)]))) ,3) ./ (sz(3)*ones([sz(1:2)])); % Two tailed
    critval = prctile(dCave_null,95,3);
            
end
    

function z = get_min_phi_diff(z)
    
    ind=z>pi; z(ind)=z(ind)-2*pi; ind=z<-pi; z(ind)=z(ind)+2*pi;        % Bounds of maximum shift should be -pi to pi
    %z = abs(z);
    
end

function [zsc] = calc_normalized_dCave(Cave1,Cave2,Cavep1,Cavep2, do_phi)
 
    % Get diff
    dCave = Cave1 - Cave2;
    if do_phi; dCave = get_min_phi_diff(dCave); end
    dCave = abs(dCave);

    % Git diff null
    dCave_null = Cavep1 - Cavep2;
    if do_phi; dCave_null = get_min_phi_diff(dCave_null); end
    sz = size(dCave_null);
    
    % Centre both values about zero (new code)
    null_means = mean(dCave_null,3);
    dCave = dCave - null_means;
    dCave_null = dCave_null - repmat(null_means,[1,1,sz(3)]);

    zsc = (dCave)./mean(abs(dCave_null),3);
    
    %pvals = sum(dCave_null > abs(repmat(dCave,[1,1,sz(3)])),3) ./ (sz(3)*ones([sz(1:2)])); % One tailed
    pvals = sum( (dCave_null > abs(repmat(dCave,[1,1,sz(3)]))) | (dCave_null < -abs(repmat(dCave,[1,1,sz(3)]))) ,3) ./ (sz(3)*ones([sz(1:2)])); % Two tailed
    critval = prctile(dCave_null,95,3);
            
end
