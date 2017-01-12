
function [data_out, xout, x2out] = extract_data(group,abscissa,do_sgolay,xld,yld,myfieldname)

    if nargin < 6
        myfieldname = 'data';
    end
    
    data = group.(myfieldname);
    xdata = group.xdata;
    if ~isempty(xdata); xout = xdata;
    else xout = abscissa;
    end
    
    x2out = group.xdata2;

    % If data still isn't real, we convert it to real, only after
    % taking mean.
    if ~isreal(data)
        data = mean(data,2);
        data = abs(data);
        %data = pls_complex2angle_centpi(data);
        warning('Complex data received! We should be bootstrapping now, meaning that input data should always be real.');
    end
        
    if ~isempty(data) && do_sgolay
        fprintf('Filtering data with SG order 3 filter.\n');
        df = mode(diff(xout));
        sgf = round(20/df);
        if mod(sgf,2) == 0; sgf=sgf+1;end
        data_out = sgolayfilt(data,3,sgf);
    else
        data_out=data;
    end
    
    % Trim trim y-axis (freq) data as needed
    %yld = group.ylims_desired;
    if ~isempty(yld)   % X data holds frequencies, which go along y-axis
        ind = xout >= yld(1) & xout < yld(2);
        xout = xout(ind);
        data_out = data_out(ind,:,:);
    end
    
    % Trim x-axis (time) data as needed
    %xld = group.xlims_desired;
    if ~isempty(xld)            % X-axis data is time, which is stored in xdata2
        ind = x2out > xld(1) & x2out < xld(2);
        x2out = x2out(ind);
        data_out = data_out(:,:,ind);
    end

end