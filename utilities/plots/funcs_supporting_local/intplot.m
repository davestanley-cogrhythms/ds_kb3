

function h = intplot (t,X)

    [Xint, Xdiff] = calc_Xint(X);
    hold off;
    scale_val = std(Xint);
    plot(t,zscale(X(:),scale_val));
    hold on; plot(t,(Xint(:)),'r','LineWidth',2);
    Xmap = 2.^(-abs(Xdiff(:)));
    hold on; plot(t,zscale(Xmap(:),scale_val)-10*scale_val,'k','LineWidth',1);
    %hold on; plot(t,zscore(Xdiff(:))-10,'k');
    hold off;

end