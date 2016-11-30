
function h=subplotsq(N,i)
%     h=subplotsq(N,i)
% 
%     Wrapper shortcut for producing subplots in arranged
%     in a square (see code below)

    Ncols = ceil(sqrt(N));
    Nrows = ceil(N/Ncols);
    h=subplot(Nrows,Ncols,i);
    
end