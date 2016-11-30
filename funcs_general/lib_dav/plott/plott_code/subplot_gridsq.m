
function varargout=subplot_gridsq(N)
%     h=subplotsq(N,i)
% 
%     Wrapper shortcut for producing subplots in arranged
%     in a square (see code below)

    Ncols = ceil(sqrt(N));
    Nrows = ceil(N/Ncols);
    [varargout{1:nargout}]=subplot_grid(Nrows,Ncols);
    
end