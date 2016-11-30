

function h = filtplot(t,X)
    plot(t,X);
    [Xout] = qif(t,X,[0 58; 62 Inf]); hold on; plot(t,Xout,'r')
    %[Xout] = qif(t,X,[58 62]); hold on; plot(t,Xout,'k','LineWidth',2)
    hold off;

end