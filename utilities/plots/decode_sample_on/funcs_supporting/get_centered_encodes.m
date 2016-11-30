

function [indg] = get_centered_encodes(ind2,Npoints)

    %Npoints = 29;
    indg = repmat(ind2(:)',Npoints,1) + repmat([[1:Npoints] - round(Npoints/2)+0]',1,length(ind2));  % points x trials



end