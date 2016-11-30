
function varargout = venn_labeled(mylabels,add_percent_sign,fixed_sizes,varargin)
% [H, S] = venn_labeled(mylabels,add_percent_sign,inputs_to_venn)
% Wrapper function for venn.m to add labels
% Varargin and varargout are the same as for venn.m (venn.m).
% Assumes venn of the form venn(A,I) and not venn(Z).
%     mylabels - either empty or a cell array of labels of length(A) for each entry in A.
%     display_percent_sign - If 1, adds a "%" to the end of each label. Empty = 0 (default)
%     fixed_sizes - If venn.m spits out errors, try setting fixed_sizes to 1. This will produce a venn diagram with default sizes.
% Dave Stanley, BostonU, stanleyd@bu.edu.


    if isempty(add_percent_sign)
        add_percent_sign = 0;
    end
    
    if isempty(fixed_sizes)
        fixed_sizes = 0;
    end
    
    A = varargin{1};
    I = varargin{2};
    
    Nsections = length(A) + length(I);
    
    if ~fixed_sizes
        [myouts{1:2}] = venn(varargin{:});
    else
        if length(A) == 3
        [myouts{1:2}] = venn([1,1,1],[1 1 1 1]*0.2);
%         maxA = min(A);
%         [myouts{1:2}] = venn(A,[maxA maxA maxA maxA]*0.2);
        elseif length(A) == 2
            [myouts{1:2}] = venn([1,1],[1]*0.2);
        end
    end
    
    S=myouts{2};
    
    AI = [A(:); I(:)]';
    
    % Add Labels and numbers
    for i = 1:length(A)
        
        mytitle = num2str(AI(i),'%5.3g');
        if ~isempty(mylabels); mytitle = [mylabels{i} ' ' mytitle]; end
        if add_percent_sign; mytitle = [mytitle '%'];end
        
        text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), mytitle,'HorizontalAlignment','center');
    end
    
    % Add just numbers
    for i = length(A)+1:Nsections
        
        mytitle = num2str(AI(i),'%5.3g');
        if add_percent_sign; mytitle = [mytitle];end
        
        text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), mytitle,'HorizontalAlignment','center');
    end

    varargout = myouts;
end


