
function z = get_circ_diff(z,is_circular)
    % Treats differences as circular measures in radians (only if
    % is_circular = 1). Otherwise, returns input unchanged
    
    if nargin < 2
        is_circular = 1;
    end
    
    if is_circular
        ind=z>pi; z(ind)=z(ind)-2*pi; ind=z<-pi; z(ind)=z(ind)+2*pi;        % Bounds of maximum shift should be -pi to pi
        %z = abs(z);
    end
    
end