
function out = get_clist(i,format)
    % Returns a list of colors for plotting
    % i specifies which color to return
    % format specifies whether to return a character or a vector (see here:
%         http://www.mathworks.com/help/matlab/ref/colorspec.html )
    
    colourarr = 'bgrcmkbgrcmkbgrcmkbgrcmkbgrcmkbgrcmkbgrcmkbgrcmkbgrcmkbgrcmk';
    colourmat = [0 0 1; ...
                 0 1 0; ...
                 1 0 0; ...
                 0 1 1; ...
                 1 0 1; ...
                 0 0 0; ...
                 ];
     colourmat = repmat(colourmat,10,1);
             
    

    if nargin < 1
        i=[];
    end
    if isempty(i)
        i=1:length(colourarr);
    end
    
     
    if nargin < 2
        format = 'strings';
    end
        
    switch format
        case 'strings'
            if isvector(i)
                out = arrayfun(@(i) colourarr(mod(i-1,length(colourarr))+1),i);
            else
                out = colourarr(mod(i-1,length(colourarr))+1);
            end
        case 'matrix'
            out = colourmat(i,:);
    
    end
end