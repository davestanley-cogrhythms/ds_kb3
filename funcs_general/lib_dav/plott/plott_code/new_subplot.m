
function [hsp, i, returns] = new_subplot(Nmax,i,hsp,use_subplot_grid, fig_handle,visible)

    if nargin < 4
        use_subplot_grid = 0;
    end
    if nargin < 5
        fig_handle = @figl;
    end

    if nargin < 6
        visible = 'on';
    end

    % i is total number of subplots in all figures
    % j  is number of subplots in current figure
    j = mod(i,Nmax);
    j(j==0) = Nmax;
    
    % Example:
    % If there are 4 subplots, then:
    %    i = [ 1     2     3     4     5     6 ]
    %    j = [ 1     2     3     4     1     2 ]
    
    
    returns = [];
    if j == 1 && i > Nmax
%         fprintf('Current figure full; resetting \n');
%         pause
        prompt = 'Type a comment or hit enter to continue q to quit: ';
        returns = input(prompt,'s');
        if strcmp(returns,'q') || strcmp(returns,'Q'); return; end
    end

    %[i j]
    if use_subplot_grid
        if j==1; fig_handle('visible',visible); hsp = subplot_gridsq(Nmax); end   % Start a new figure if necessary
        hsp.set_gca(j);
    else
        if j==1; fig_handle();end                                         % Start a new figure if necessary
        subplotsq(Nmax,j);
    end
    
    
    
    i=i+1;


end