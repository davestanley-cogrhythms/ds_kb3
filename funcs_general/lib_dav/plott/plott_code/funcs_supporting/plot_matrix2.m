

function handles = plot_matrix2 (t, data, opt_strct, varargin)



    ds = 1;
    shift = 0;
    zero_means = 0;
    plotloglog = 0;
    normalize_everything = 0;

    if exist('opt_strct','var')
        if ~isempty(opt_strct);
            if isfield (opt_strct,'ds'); ds = opt_strct.ds; end
            if isfield (opt_strct,'shift'); shift = opt_strct.shift; end
            if isfield (opt_strct,'zero_means'); zero_means = opt_strct.zero_means; end
            if isfield (opt_strct,'plotloglog'); plotloglog = opt_strct.plotloglog; end
            if isfield (opt_strct,'normalize_everything'); normalize_everything = opt_strct.normalize_everything; end
        end
    end

    if normalize_everything
       for jj = 1:size(data,2)
           xtemp = data(:,jj);
           xtemp = xtemp - mean(xtemp);
           xtemp = xtemp / std(xtemp);
           data(:,jj) = xtemp;
       end
    end

    if ds > 1
        t = downsample(t, ds);
        data = downsample(data, ds);
    end

    for i = 1:size(data,2)
        data(:,i) = data(:,i)+shift*(i-1);
    end
    handles = plot(t,data,varargin{:}); 



end