


function [hl1,hl2] = plott_matrix3D(varargin)
%     [hl,hl2] = plott_matrix3D(t,X,'options')
%     Flexible plotting for 3D matrices
%     Uses the boundedline function (by Kelly Kearney, available from Mathworks)
%     FORMS
%         [hl1,hl2] = plott_matrix3D(X)
%         [hl1,hl2] = plott_matrix3D(t,X)
%         [hl1,hl2] = plott_matrix3D(X,'options')
%         [hl1,hl2] = plott_matrix3D(t,X,'options')
%     INTPUTS
%         X - vector or matrix of data. Can be 2D or 3D. If 3D, then one of
%         the dimensions will be used for averaging. 
%         options - specifies options in the form of name and value pairs
%             (i.e. 'fs',1024)
%     OPTIONS (name value pairs)
%         'randcol', boolean -  randomize columns for plotting. BOOLEAN = 1 for true
%         'fname', fname - function handle or cell array of multiple function handles
%         'fs' - sampling freq of data
%         'active_dim' - sets the active dimension along which operations (eg do_mean, do_shift, etc)
%               are performed. If both do_mean and do_shift are active, do_mean will take priority
%               and act on the active dimension, while do_shift will operate on the next highest
%               dimension.
%         'do_mean' - if 1, will take the mean along active_dim
%         'showErrorbars' - if do_mean=1, will plot error bars (mean +/- SEM)
%         'do_shift' - Will shift entries of active_dim by this value, effectively spreading out the data
%         'do_zscore' - applies z-score to 1st dimension of data
%         'do_filter' - A vector of length size(X,active_dim). Performs logical indexing on active_dim
%             Useful for removing bad data.
%         'Color' - Color arguments of the form accepted by plot LineSpec (i.e. 'b-','kx',etc.)
%         'LineSpec' - takes in a cell array of linespec arguments, passes these to the plot command
%         'LineSpec2' - takes in a cell array of linespec arguments, passes
%             these to the boundedline command.
%     OUTPUTS
%         [hl1] - plot handle of line
%         [hl2] - plot handle of bounded line
% 
%     EXAMPLES
%         t = linspace(0, 2*pi, 50);
%         X1 = sin(t);
%         X2 = 0.5*cos(t);
%         
%         X1 = repmat(X1(:),1,100);
%         X2 = repmat(X2(:),1,100);
%         E1 = randn(size(X1,1),100)*.5;
%         E2 = randn(size(X2,1),100)*0.2;
% 
%         X = cat(3, X1+E1, X2+E2);
% 
%         figure;
%         subplotrows(2,1); plot(X(:,:,1)); title('Group 1 - basic plot');
%         subplotrows(2,2); plot(X(:,:,2)); title('Group 2 - basic plot');
% 
%         figure; title('Averaged with error bars showing SEM (standard error in mean)');
%         h=plott_matrix3D(t, X,'do_mean',1,'do_zscore',0,'showErrorbars',1,'LineSpec',{'-'})
% 
%         figure; title('All data, without mean. Grouped along 2nd dimension');
%         h=plott_matrix3D(t, X,'do_mean',0,'do_zscore',0,'showErrorbars',0,'LineSpec',{':'})
% 
%         figure; title('Averaged with error bars showing SEM (standard error in mean)');
%         h=plott_matrix3D(t, X,'do_mean',0,'do_zscore',0,'showErrorbars',0,'do_shift',10,'LineSpec',{'-'})
% 
%         figure; title('Averaged with error bars showing SEM (standard error in mean)');
%         h=plott_matrix3D(t, X,'do_mean',0,'do_zscore',0,'showErrorbars',0,'do_shift',10,'active_dim',3,'LineSpec',{'-'})
%         
%         myfilter = [false(1,97) true(1,3)];        % Fitler to use only 3 of the original 100 data points in our average
%                                                    % Useful for removing bad data
%         figure; title('Averaged with error bars showing SEM (standard error in mean)');
%         h=plott_matrix3D(t, X,'do_mean',1,'do_zscore',0,'showErrorbars',1,'LineSpec',{'-'},'do_filter',myfilter)
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
%         March 27 2014 - Updated to allow t to be a vector (Dave)



    
    [reg, args]=parseparams(varargin);  % Regs should contain t and X    
    p = inputParser;
    addOptional(p,'fs',1,@isnumeric);
    addOptional(p,'showErrorbars',1,@isnumeric);
    addOptional(p,'active_dim',2,@isnumeric);
    addOptional(p,'do_mean',1,@isnumeric); % Means along active dim
    addOptional(p,'do_shift',[],@isnumeric);    % Shifts along active dim
    addOptional(p,'do_filter',[],@islogical);    % Filters along active dim
    addOptional(p,'randcol',0,@isnumeric);    % Filters along active dim
    addOptional(p,'do_zscore',0,@isnumeric);
    addOptional(p,'Color',{},@iscell);
    addOptional(p,'LineSpec',{},@iscell);
    addOptional(p,'LineSpec2',{},@iscell);
    addOptional(p,'fname',@plot);
    
    parse(p,args{:});
    vars_pull(p.Results)
    
    % Set Defaults
    if isempty(fs); fs = 1; end
    
    [t,X,fs] = sort_out_timeseries(reg,fs); % Assigns values to t, X, and fs. If available, the fs from t always overrides fs.

    
    X = activedim_shift(X,active_dim);  % Shift matrix dimensions so dim 2 is always active
    
    if ~isempty(do_filter)
        X=X(:,do_filter,:);
    end
    
    if do_zscore
        X=zscore(X);
    end
    
    if size(X,2) <= 1   % If only one dimension, there is nothing to mean!
        showErrorbars = 0;
        do_mean = 0 ;
    end
    
    if showErrorbars; do_mean=1; end    % If showing error bars, need to take mean
    
    if do_mean
        Xste=std(X,[],2) / sqrt(size(X,2)); % Standard error
        %Xste=std(X,[],2);                   % Standard deviation
        X=mean(X,2);
        %X=median(X,2);
        t = mean(t,2);      % Some extra operations, but it's worth it.
    end
    
    if ~isempty(do_shift)
        permuted = 0;
        if size(X,2) <= 1; X = permute(X,[1,3,2]); permuted = 1; end  % If already meaned dim2, bump up next dimension
        shiftmat = (0:size(X,2)-1)*do_shift;
        shiftmat = repmat(shiftmat,[size(X,1),1,size(X,3)]);
        X=X+shiftmat;
        if permuted; X = permute(X,[1,3,2]);end  % Return to defualt
    end

    N = size(X,1);
    Nepochs = size(X,2);
    
    if randcol cols = randperm(Nepochs);
    else cols = 1:Nepochs; end
    
    for i = 1:length(cols)
        hl1{i} = fname(t(:,cols(i)),squeeze(X(:,cols(i),:)),Color{:},LineSpec{:}); hold on; 
        
        hl2 = [];
        if showErrorbars
            
            % Transfer cmap over to LineSpec2 for boundedline if not already present
            if isempty(find(strcmp(LineSpec2,'cmap')))
                lsc_ind = find(strcmp(LineSpec,'Color'));
                if ~isempty(lsc_ind)
                    if ~ischar(LineSpec{lsc_ind+1})                 % Linespec is of form 'Color',[0 1 0]
                        LineSpec2{end+1}='cmap'; % The color specification is a row vector
                        LineSpec2{end+2}=LineSpec{lsc_ind+1};
                    else
                        Color=LineSpec(lsc_ind+1);         % Linespec is of form 'Color','b-'
                    end
                end
            end

            
            
            hold on; hl2{i} = boundedline(t,squeeze(X),permute(squeeze(Xste),[1 3 2]),Color{:},'alpha',LineSpec2{:});
            
            %hold on; hl{i} = shadedErrorBar(repmat(t(:),1,size(squeeze(X),2)),squeeze(X),squeeze(Xste));
            %outlinebounds (hl{i});
        end
    end
    hold off;
    
    if length(hl1) == 1
        %hl1 = cell2mat(hl1);
        hl1 = hl1{1};
    end
    if length(hl2) == 1
        %hl2 = cell2mat(hl2);
        hl2 = hl2{1};
    end
    
%     if do_zscore; ylim([-1.5 2])
%     end
    
end


function vars_pull(s)
    for n = fieldnames(s)'
        name = n{1};
        value = s.(name);
        assignin('caller',name,value);
    end
end


