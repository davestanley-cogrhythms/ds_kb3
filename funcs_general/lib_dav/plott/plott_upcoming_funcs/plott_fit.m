

function [R, P, rsq2] = plott_fit(varargin)
    % Plot and fit an order 1 polynomial

    use_normalization = 1;

    x = varargin{1};
    y = varargin{2};
    x=x(:);
    y=y(:);
    
    if any(size(x) ~= size(y)); error('Vectors must be the same length.'); end
    
    [~,ind] = sort(x);
    x=x(ind);
    y=y(ind);
%     
%     x=zscore(x);
%     y=zscore(y);
%     
    htemp = plot(x,y,varargin{3:end});
    currcolour = get(htemp,'Color');
    
    [R, P] = corrcoef(x,y);
    if use_normalization
        [p,S,mu] = polyfit(x,y,1);
        [y0,delta0] = polyval(p,x,S,mu);
    else
        [p,S] = polyfit(x,y,1);
        [y0,delta0] = polyval(p,x,S);
    end
    
    

    yresid = y - y0;
    SSresid = sum(yresid.^2);
    SStotal = (length(y)-1)*var(y);
    rsq = 1 - SSresid/SStotal;
    
    normr = S.normr;
    rsq2 = 1-normr^2/SStotal;
    
    
    hold on; boundedline(x,y0,delta0,'alpha','cmap',currcolour)
    %hold on; errorbar(x,y0,delta0)
    %[R, P] = corrcoef(mydata,randn(size(mydata)))
    
    title(['Corr coef ' num2str(R(2)) ' p=' num2str(P(2)) ' N=' num2str(length(x)) ' slope=' num2str(p(1))],'FontSize',22);
    
    %hold on; hl{i} = boundedline(t,squeeze(X),permute(squeeze(Xste),[1 3 2]),'alpha',LineSpec{:});
    
    legend('data','50% error');
    
end