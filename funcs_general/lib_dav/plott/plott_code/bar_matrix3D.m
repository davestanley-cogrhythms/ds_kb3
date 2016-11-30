

function [hbar herr] = bar_matrix3D(varargin)
%     h = bar_matrix3D(varargin)
%     Takes in a 3D matrix of data, calculates means, standard deviations,
%     and significances across the columns (by default) and returns a bar
%     plot with significance indicators
%     Uses function errbar_sig. Also uses  the plotwitherr function (available 
%     from Mathworks, by Martina F. Callaghan)
%     FORMS
%         [hbar herr] = bar_matrix3D(X)
%         [hbar herr] = bar_matrix3D(n,X)
%         [hbar herr] = bar_matrix3D(X,options,bar_vars)
%         [hbar herr] = bar_matrix3D(n,X,options,bar_vars)
%     INTPUTS
%         n - coordinates along abscissa
%         X - vector or matrix of data
%         options - specifies options in the form of name and value pairs
%         bar_vars - variables as passed to conventional bar plot command
%     OPTIONS
%         'active_dim' - dimension along which to conduct statistics
%         'statfun' - handle for a function that returns a hypothesis (h= 1 or 0).
%             This will be used to conduct statistics relative to 1st datapoint.
%         'mysymbol' - symbol to use for denoting significance
%         'myfontsize' - Font size of symbol (default 25)
%         'star_shift_factor' - amount by which to shift the significance "star" 
%     OUTPUTS
%         [hbar herr] - plot handles
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 
    
    [reg, args]=parseparams(varargin);  % Reg should contain X and maybe n
    [statfun, args] = parse_pair(args,'statfun',@myttest);
    [active_dim, args_leftover] = parse_pair(args,'active_dim',[]);
    
    if length(reg) == 1
        X = reg{1};
        n = 1:size(X,1);
    elseif length(reg) == 2
        X = reg{2};
        n = reg{1};
    else
        fprintf('Incorrect number of input arguments \n');
        return
    end
    
    Xdims = ndims(X);
    if isempty(active_dim)
        if Xdims == 2
            active_dim = 1;
        elseif Xdims == 3
            active_dim = 2;
        else
            fprintf('X must be 2D or 3D. Exiting \n'); return;
        end
    end
    
    X = activedim_shift(X,active_dim);  % Shift matrix dimensions so dim 2 is always active
    
    Xmu = squeeze(mean(X,2));
    Xste=squeeze(std(X,[],2) / sqrt(size(X,2)));
    
    sz = size(X);
    
    if Xdims == 2
        h = false(sz(1),1);
        p = zeros(sz(1),1);
        for i = 2:sz(1)
            h(i) = statfun(X(i,:),X(1,:)); % Stats along COLS
            p(i) = myttestp(X(i,:),X(1,:)); % p value is not used...for now
        end
    elseif Xdims == 3
        h = false(sz(1),sz(3));
        p = zeros(sz(1),sz(3));
        for i = 1:sz(1)
            for j = 2:sz(3)
                h(i,j) = statfun(X(i,:,j),X(i,:,1));
                p(i,j) = myttestp(X(i,:,j),X(i,:,1)); % p value is not used...for now
            end
        end
    else
        printf('Wrondg dimensionality of X \n');
        return;
    end
    
    [hbar herr] = bar_errsig(h,Xste,Xmu,args_leftover{:});
    
end

function h = myttest(varargin)
    %alpha = 0.05;
    %[h ~] = ttest(varargin{:},alpha);
    [h, ~] = ttest(varargin{:});
    %[~, h] = signrank(varargin{:});

end

function p = myttestp(varargin)
    %alpha = 0.05;
    %[h p] = ttest(varargin{:},alpha);
    [~, p] = ttest(varargin{:});
    %[p, h] = signrank(varargin{:});
end


