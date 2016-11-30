

function y = replicate_inline(x,Narr)
    % Takes in a vector x of length N and a vector Narr, also of length N
    % For each entry x(i), this function drops Narr(i) copies of x(i) in
    % the output y. It then moves on to the next entry in x.
    % 
    % For example x=[5 2 6] and Narr = [3 2 1] would produce
    % y = [5 5 5 2 2 6]
    %
    % In the context of my code...
    % This is used for taking arrays that are describing values associated
    % with files and expanding them to have a single entry for every cell
    % described by that file.
    
    if ~isvector(x); warning('x must be vector'); end
    if ~isvector(Narr); warning('Narr must be vector'); end
    %ic = iscolumn(x);
    
    x = x(:); Narr = Narr(:);
    x = num2cell(x);
    Narr = num2cell(Narr);
    y = cellfunu(@(x,y) repmat(x,y,1),x,Narr);
    y = vertcat(y{:});
    

end