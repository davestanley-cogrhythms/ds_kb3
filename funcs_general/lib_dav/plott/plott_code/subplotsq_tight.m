
function h=subplotsq_tight(N,i)
%     h=subplotsq(N,i)
%     requires tight_subplot.m from Mathworks central
% 
%     Wrapper shortcut for producing subplots in arranged
%     in a square (see code below)

    Ncols = ceil(sqrt(N));
    Nrows = ceil(N/Ncols);
    %gap_default_val = [0.01,0.01];
    gap_default_val = [0.06,0.06];
    h=subtightplot(Nrows,Ncols,i,1*gap_default_val);
    
end