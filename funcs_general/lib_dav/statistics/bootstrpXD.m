

function [BOOTSTAT, BOOTSAM] = bootstrpXD(varargin)
%     BOOTSTAT = bootstrpXD(dim,NBOOT,BOOTFUN,D1,...)
%     Does bootstrapping with X dimensional inputs
% 
%     OVERVIEW
%     With dim=[], bootstrpXD operates exactly the same the familiar Matlab
%     bootstrp function, except the matricies D1,... can be greater than
%     2D. 
%     Assigning a scalar value to dim tells tells bootstrpXD which dimension to
%     shuffle along ([] for default, dim=1, shuffle rows).
% 
%     The formats are the same as bootstrp, except for the first argument,
%     dim. See bootstrp documentation.
%     BOOTSTAT = bootstrp(dim,NBOOT,BOOTFUN,D1,...)
%     [BOOTSTAT,BOOTSAM] = bootstrp(dim,...)
%     BOOTSTAT = bootstrp(dim,...,'Name',Value) 
%     BOOTSTAT = bootstrp(dim,..., 'PARAM1',val1, 'PARAM2',val2, ...) specifies
%         optional parameter name/value pairs to control how bootstrp performs
%         computations.  Parameter names/values may only appear after the data
%         arguments used as inputs to BOOTFUN.  Parameters are:
% 
%            'Weights' - Observation weights
% 
%            'Options' - options including 
%                         'UseParallel'
%                         'UseSubstreams'
%                         'Streams'

%     IMPORTANT NOTE
%     The return variable of the function handle BOOTFUN must be a scalar
%     or a vector. Matrix outputs will be reshaped to column vectors.
% 
%     ALGORITHM
%     if D1 is ND dimensional, it will first permute D1 to place dim, the
%     shuffling dimension of interest, first (rows). Then it will pack the
%     remaining dimensions into the columns, and call Matlab's bootstrp command
%     using a wrapper function. This wrapper function will undo all of the
%     aforementioned transformations prior to calling the user's input
%     function, BOOTFUN. See documentation for bootstrp for more details.
% 
% 
%     % EXAMLES
%     
%     %% Simpler example
%     % Generate some data
%     x = 1:(55*7*2);
%     x = reshape(x,[55,7,2]);
%     x = x + 0.1*mean(x(:))*randn([55,7,2]) - mean(x(:));
%     
%     % Define some function
%     myfunc = @(x) sum(x,2); % Sum along columns
%     xout = myfunc(x);
%     
%     % Plot original data
%     figure; h1 = plot(squeeze(xout),'k','LineWidth',4);
%     
%     % Plot bootstrapped data
%     y = bootstrpXD(2,300,myfunc,x);     % Shuffle data along dim 2
%     y = reshape(y,[300,55,2]);   % Unpack bootstrapped data; shuffles along 1st dimension
%     ybar = mean(y);
%     ystd = std(y);
%     hold on; h2 = errorbar(squeeze(ybar),squeeze(ystd),'r.','MarkerSize',20,'LineWidth',1)
%     legend([h1(1),h2(1)],'Data','Bootstrap Mean & Standard Deviation')
%     
%     
%     %% Higer dimensional example
%     % Generate some high dimensional input
%     x = 1:(6*2*3*4*5);
%     x = reshape(x,[6,2,3,4,5]);
%     x = x + 0.1*mean(x(:))*randn([6,2,3,4,5]) - mean(x(:));
%     sz = size(x);
%     
%     % Define some function
%     myfunc = @(x) sum(x);
%     xout = myfunc(x);
%     
%     % Plot original data
%     figure; h1=plot(squeeze(xout(:,2,:,:,5)),'k','LineWidth',4);
%     
%     % Plot bootstrapped data
%     y = bootstrpXD([],30,myfunc,x);
%     y = reshape(y,[30,2,3,4,5]); % Unpack y
%     ybar = mean(y);
%     ystd = std(y);
%     hold on; h2=errorbar(squeeze(ybar(:,2,:,:,5)),squeeze(ystd(:,2,:,:,5)),'r.','MarkerSize',20,'LineWidth',1)
%     legend([h1(1),h2(1)],'Data','Bootstrap Mean & Standard Deviation')
%     
%     % David Stanley, Boston University

    [reg, prop]=parseparams(varargin);
    i=0;
    i=i+1; dim=reg{i};
    i=i+1; N=reg{i};
    i=i+1; myfunc=reg{i};
    i=i+1; mydata=reg(i:end);
    
    if isempty(dim); dim=1;end
    
    for i = 1:length(mydata)
        x = mydata{i};
        if ~isscalar(x)
            sz = size(x);
            Nd = ndims(x); x = permute(x,[dim, 1:dim-1 dim+1:Nd]);  % Brings dimension of interest to the front (i.e. rows)
                         % To test this code, run>> dim=2; x = randn([1 2 3 4 5 6]); Nd = ndims(x); x = permute(x,[dim, 1:dim-1 dim+1:Nd]); size(x)
            szp = size(x); cpsz = cumprod(szp); szrest = cpsz(end) ./ szp(1);      % Total amount of data in remaining dimensions
            x = reshape(x,[szp(1),szrest]);
            mydata{i} = x;      % Pack data back into cell array
            
            % Pack useful metadata (necessary to undo the transform)
            s.sz{i} = sz;     
            s.szp{i} = szp;
            s.Nd{i} = Nd;
            s.dim = dim;
            
        end
    end
    
    s.myfunc = myfunc;  % Also pack function handle into structure.

    % Call Matlab's bootstrp function with reshaped data. The wrapper
    % function, myfunc2 will undo the above transforms prior to calling the
    % input function BOOTFUN.
    [BOOTSTAT, BOOTSAM] = bootstrp(N,@myfunc2,s,mydata{:},prop{:});


end


function out = myfunc2(varargin)
% Undo all the transforms before running myfunc

    % Unpack inputs
    s = varargin{1};            % Input metadata structure
        sz = s.sz;
        szp = s.szp;
        myfunc = s.myfunc;
        Nd = s.Nd;
        dim = s.dim;
    mydata = varargin(2:end);   % Data (which will be shuffled)
    
    % Undo all the transformations
    for i = 1:length(mydata)
        x = mydata{i};
        if ~isscalar(x)
            % Undo reshape
            x = reshape(x,szp{i});
            
            % Undo permutation
            permat = 1:Nd{i};
            permat([1,dim]) = [dim,1];  % Swap dim and 1
            x = permute(x,[2:dim 1 dim+1:Nd{i}]);      % Put dim back in place
%                 To test this code run>>
%                     dim=5; x = randn([1 2 3 4 5 6]); Nd = ndims(x); x = permute(x,[dim, 1:dim-1 dim+1:Nd]); size(x)    % Do permute
%                     x = permute(x,[2:dim 1 dim+1:Nd]); size(x)
            
            % Repack data
            mydata{i} = x;
        end
    end
    

    % Run the function
    out = myfunc(mydata{:});
    
    % Make sure output is a row
    if isvector(out)
        out = out(:)';
    end

end


