


function [returns] = plott_ani_pairs(varargin)
%     [returns] = plott_ani_pairs(X1,fhandle1,X2,fhandle2,...XN,fhandleN,options)
%     Takes a set of matrices X1,X2,X3, etc. paired with an equal number
%     of function handles. Creates a figure with a bunch of subplots, and runs
%     plott_ani for each of the subplots across the columns of X1,X2,X3. See
%     examples for clarity.
%     FORMS
%         [returns] = plott_ani_pairs(X1,fhandle1,X2,fhandle2,...XN,fhandleN,options)
% 
%     INTPUTS
%         Xi - vector or matrix of data, where rows are times and columns are variables
%            - if Xi is 3D, then the 2nd dimension is used for animation and the 3rd is plotted with hold on.
%         fhandlei - function handle to be paired with Xi
%         options - specifies some plotting options in the form of name and value pairs
%     OPTIONS
%         'randcol', boolean -  randomize columns for plotting. BOOLEAN = 1 for true
%         'splitup',N - automatically converts large datafiles into smaller chuncks of size N that are
%             more agreeable for plotting. Zero for no splitting (default).
%         'plotargs' - cell array of option pairs to be passed to all plotting functions
%         'fsubplot' - subplot function to use. Default is @subplotsq
% 
%     OUTPUTS
%         returns - set of user keystrokes
% 
%     EXAMPLE
%     Example 1 - Plot 3 separate 3D matrices
%     X = abs(randn([100,5,2])); Y = -X; Z = (X).^2;
%     out = plott_ani_pairs(X,@plot,Y,@plot,Z,@loglog,'fsubplot',@subplotsq);
%     
%     Example 2 - Plot 3 separate 2D matrix that are converted to 3D using the 
%                 splitup option
%     fs = 0.9;
%     t=(1:1/fs:950)'; X = [sin(2*pi/21*t) 1.1*sin(2*pi/21*t)+0]; Y = -X; Z = abs(X)+1;
%     out = plott_ani_pairs(X,@plot, ...
%         Y,@(x) plot([1:length(x)]/fs,x), ...
%         Y,@(x) plott_matrix3D(x,'fs',fs,'do_mean',0,'do_shift',0.3), ...
%         Z,@loglog,'fsubplot',@subplotrows,'splitup',100);
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 

    [reg, args]=parseparams(varargin);  % Regs should contain t and X
    p = inputParser;
    addOptional(p,'randcol',0,@isnumeric);
    addOptional(p,'splitup',0,@isnumeric);
    addOptional(p,'fsubplot',@subplotrows);
    addOptional(p,'plotargs',{});
    
    parse(p,args{:});
    vars_pull(p.Results)
    
    Nreg = length(reg);
    reg = reshape(reg,2,round(Nreg/2))';
    fhandles = reg(:,2);
    data = reg(:,1);
    
    Nplots = size(reg,1);
    
    if ~dimcheck(data); fprintf('# columns in inputs must all be the same. \n'); return; end
    
    X = data{1}; % Take prototypical data pair.
    
    % Split up data for the heck of it if the input is a column
    if isvector(X)
        for i = 1:Nplots
            data{i}=data{i}(:);
            if splitup <= 0 splitup = round(length(X)/10); end % Turn on splitup, since it's pointless otherwise
        end
    end
    
    if splitup > 0
        for i = 1:Nplots
            data{i} = splitup_data(data{i},splitup);     % Take the data in each row, and split it into columns, which will then be animated across
            if isnan(sum(data{i}(:))); returns=NaN;return; end % splitup returns NaN if X is wrong dimensionality.
        end
    end

    
    Nepochs = size(data{1},2);
    
    if randcol cols = randperm(Nepochs);
    else cols = 1:Nepochs; end
    
%     for i = 2
    for i = 1:Nepochs
        for j = 1:Nplots
            fsubplot(Nplots,j);
            fhandles{j}(squeeze(data{j}(:,cols(i),:)),plotargs{:});
        end
        
        prompt = 'Type a comment or hit enter to continue q to quit: ';
        returns{i} = input(prompt,'s');
%         returns{i}=1;

        if strcmp(returns{i},'q') || strcmp(returns{i},'Q'); break; end
    end
end


function vars_pull(s)
    for n = fieldnames(s)'
        name = n{1};
        value = s.(name);
        assignin('caller',name,value);
    end
end


function are_same = dimcheck (data)
    % Checks dimensions of X to make sure they're all consistent
    
    out = cellfun(@(x) size(x,2),data);
    are_same = ~sum(~(mode(out) == out)); % Make sure that number of columns in each input are equal to the mode of all inputs

end



