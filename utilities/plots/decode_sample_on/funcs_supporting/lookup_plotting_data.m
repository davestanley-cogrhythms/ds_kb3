

function [xsq, tcent] = lookup_plotting_data(indg,ind2,x,t)
    Npoints=size(indg,1);
    xsq = x(indg);
    tsq = t(indg);
    tcent = tsq - repmat(t(ind2)',Npoints,1);
end

