

function R = plot_corrfreq(X,fs,interval)
    % R = plot_corrfreq(X,fs,interval)
    % Filters the data in the specified frequency interval, then takes the
    % correlation coefficient. Plots it if plot_on == 1
    
    plot_on=1;

    if isvector(X);X=X(:);end
        
    Nrows = size(X,1);
    Ncols = size(X,2);
    t = 1:Nrows; t=t/fs;
    
    Xfilt = zeros(Nrows,Ncols);
    
    for i = 1:Ncols
        Xfilt(:,i) = qif(t,X(:,i),interval);
    end
    
    R = corrcoef(Xfilt);
    if plot_on   
        imagesc(R);
        colorbar
    end

end



%Quick ideal filter
function dout = qif (datatimes, data, interval)

    ts1 = timeseries(data,datatimes);
    ts1filt = idealfilter (ts1, interval, 'notch');
    dout = ts1filt.data;
    
end