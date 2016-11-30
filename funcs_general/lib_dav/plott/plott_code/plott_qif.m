

function hl = plott_fs(varargin)
%     hl = plott_fs(varargin)
%     Like normal plot, but takes in an argument 'fs' and uses this to
%     generate a time vector. This is designed to be compatible with the
%     format used by other plotting commands in this package.
%     FORMS
%         [hl] = plott_fs(x,...)            
%         [hl] = plott_fs(t,x,...)          
%         [hl] = plott_fs(x,'fs',fs,...)    
%         [hl] = plott_fs(t,x,'fs',fs,...)  
%     INTPUTS
%         t - input times, passed on to conventional plot function (spacing of t takes
%             priority over fs if both are passed)
%         x - input data, passed on to conventional plot function
%         fs - sampling rate (Default=1; when t is supplied fs is overridden by spacing of t )
%     OPTIONS
%         ... - name and value pairs that are passsed on to conventional plot
%     OUTPUT
%         hl - set of user keystrokes
    
    [reg, args]=parseparams(varargin);  % Regs should contain X and maybe t 
    [fs, args] = parse_pair(args,'fs',1);
    [myinterval, args] = parse_pair(args,'interval',[]);
    [mycolors, args] = parse_pair(args,'mycolors','b');
    [axis_lims, args] = parse_pair(args,'axis_lims',[NaN NaN]); % Remove any axis_lims arguments they may receive
    [t,X,fs] = sort_out_timeseries(reg,fs); % Assigns values to t, X, and fs. If available, the fs from t always overrides fs.
    
    Xfilt = qif(t,X,myinterval);
    
    hl = plot(t,X,mycolors,args{:});
    hold on;
    plot(t,Xfilt,'k','LineWidth',1,args{:});
    
end


%Quick ideal filter

function dout = qif (datatimes, data, interval)

    ts1 = timeseries(data,datatimes);
    ts1filt = idealfilter (ts1, interval, 'notch');
    dout = ts1filt.data;
    
end