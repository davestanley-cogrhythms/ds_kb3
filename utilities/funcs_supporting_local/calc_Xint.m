
function [Xint, Xdiff] = calc_Xint(Xvect,plot_on)

    if ~exist('plot_on','var') plot_on=0; end
    dt=get_dt;

    time_wind = 0.100; % in ms
    Xfilt = filt_amplify_linearities(Xvect);      % Diff followed by transform X -> 2^-abs(X)
    Xint = integrate_sig(Xfilt(:),round(time_wind/dt),dt);
    if nargout > 1
        Xdiff = diff_pad(Xvect,2,1,0);
        Xdiff = Xdiff(:);
    else
        Xdiff=0;
    end
    
    if plot_on
        figure; h0 = plot(zscore(Xvect(:)));
        hold on ; h1 = plot(zscore(Xint(:)),'r','LineWidth',2)
        hold on ; h2 = plot(zscore(Xfilt(:))-10,'k','LineWidth',1)
        %hold on ; h3 = plot(zscore(Xdiff(:))+10,'g','LineWidth',1)
        %legend([h0,h1,h2,h3],{'Sig','int','filtered','diff'});
        legend([h0,h1,h2],{'Sig','int','filtered'});
        figure; hist(Xint,100)
    end
end





function Y=filt_amplify_linearities(X)

    order = 2;
    X=diff_pad(X,order,1,Inf); % Pad ends with infinities, so that they will transform to zeros
    
% %     X=X+1*randn(size(X));
% %     X=1./X;
    
    %Y=X;
    Y=(2).^-abs(X);


    %Y=Y==0;
    
%     if sum(isinf(Y(:))) > 0
%         Y(isinf(Y))=max(Y(~isinf(Y)));
%         fprintf('Inf values stil found');
%     end

end



function xint = integrate_sig(x,filtsize,dt)

    if (~exist('dt','var')); dt=1; end
    
    N = length(x);
%     x = cumsum(x);
    x = abs(x);
    
    %filterwind = ones(filtsize,1)/filtsize;
    filterwind = gaussf(1:filtsize,[filtsize/2 filtsize/4]);
    
    
    xfilt = conv(x,filterwind);
    xfilt = wkeep(xfilt,N,'c');
    
%     t = (0:length(x))*dt;
%     xfilt = qif(t,x,[100 Inf]);
    
    xint = xfilt - min(xfilt);
    
%     figure;plot((0:filtsize-1)*dt,filterwind)
    

end

function y = gaussf(x,params)

    mu = params(1);
    sig = params(2);

    y = 1/sig/sqrt(2*pi)*exp(-((x-mu)/sig).^2 *1/2);

end
