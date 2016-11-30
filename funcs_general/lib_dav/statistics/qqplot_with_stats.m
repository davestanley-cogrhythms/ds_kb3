

function h = qqplot_with_stats(varargin)

    x=varargin{1};
    h = qqplot(varargin{:});

    [h, p, ksstat, cv] = kstest(x);
    warning('off','stats:adtest:OutOfRangePLow');   % Turn off warnings temporarily
    [h2, p2, adstat, cv2] = adtest(x);
    warning('on','stats:adtest:OutOfRangePLow')
    title(['N= ' num2str(length(x)) ]);
    legend( [' KS Stat=' num2str(ksstat) ' p=' num2str(p)], [ ' AD Stat=' num2str(adstat) ' p=' num2str(p2) ])

end