

function group = get_group_overrides(fname_suffix, N_criteria, gr_submode)
    % See get_group_overrides2 for simplified version; used in newer code.

    if ~exist('gr_submode','var'); gr_submode = []; end
    
    group = []; % Default output

    %N_criteria = 3;
    gr_all = Grp;
    gr_all.criteria = [2*ones(1,N_criteria)]; gr_all.criteria_alt = [2*ones(1,N_criteria)]; gr_all.ctgs = 1;


    % Mode 2
    % % Groups for reaction times; overrides others
    if ~isempty(strfind(fname_suffix,'RT')) && isempty(strfind(fname_suffix,'RT2'))
        clear group
        i=0;
        i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = 'Fast RT';
        i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = 'Slow RT';            

    end

    % Mode 3
    if ~isempty(strfind(fname_suffix,'RT2'))
        clear group
        i=0;
        if isempty(gr_submode); gr_submode = 1; end
        if any(gr_submode == 1)
            i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = 'Congruent Fast RT';
                i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = 'Congruent Slow RT';
        end
        if any(gr_submode == 2)
            i=i+1; group(i) = gr_all; group(i).ctgs = [3]; group(i).legend = 'Incongruent Fast RT';
                i=i+1; group(i) = gr_all; group(i).ctgs = [4]; group(i).legend = 'Incongruent Slow RT';
        end
           
    end

    % Mode 1
    if ~isempty(strfind(fname_suffix,'sw'))
        i=0; clear group
        i=i+1; group(i) = gr_all; group(i).ctgs = [5]; group(i).legend = 'Non-switch Trials';
        i=i+1; group(i) = gr_all; group(i).ctgs = [6]; group(i).legend = 'Switch Trials';
    end

    % Mode 4
    if ~isempty(strfind(fname_suffix,'bndry'))
        i=0; clear group
        if isempty(gr_submode); gr_submode = 1; end
        if any(gr_submode == 1)
            i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = '60% rel.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = '100% rel.';
        end
        if any(gr_submode == 2)
            i=i+1; group(i) = gr_all; group(i).ctgs = [4]; group(i).legend = '60% irr.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [3]; group(i).legend = '100% irr.';
        end
        if any(gr_submode == 3)
            i=i+1; group(i) = gr_all; group(i).ctgs = [6]; group(i).legend = '50% rel.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [5]; group(i).legend = '100% rel.';
        end
        if any(gr_submode == 4)
            i=i+1; group(i) = gr_all; group(i).ctgs = [8]; group(i).legend = '50% irr.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [7]; group(i).legend = '100% irr.';
        end
        if any(gr_submode == 5)
            i=i+1; group(i) = gr_all; group(i).ctgs = [10]; group(i).legend = '50 or 60% rel.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [9]; group(i).legend = '100% rel.';
        end
        if any(gr_submode == 6)
            i=i+1; group(i) = gr_all; group(i).ctgs = [12]; group(i).legend = '50 or 60% irr.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [11]; group(i).legend = '100% irr.';
        end
    end
    
    % Mode 5
    if ~isempty(strfind(fname_suffix,'bndABrel'))
        i=0; clear group
        if isempty(gr_submode); gr_submode = 1; end
        
        
        if any(gr_submode == 1)
            i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = '50% SchA';
            i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = '100% SchA';
        end
        if any(gr_submode == 2)
            i=i+1; group(i) = gr_all; group(i).ctgs = [4]; group(i).legend = '50% SchB';
            i=i+1; group(i) = gr_all; group(i).ctgs = [3]; group(i).legend = '100% SchB';
        end
        if any(gr_submode == 3)
            i=i+1; group(i) = gr_all; group(i).ctgs = [6]; group(i).legend = '50% Irr SchA';
            i=i+1; group(i) = gr_all; group(i).ctgs = [5]; group(i).legend = '100% Irr SchA';
        end
        if any(gr_submode == 4)
            i=i+1; group(i) = gr_all; group(i).ctgs = [8]; group(i).legend = '50% Irr SchB';
            i=i+1; group(i) = gr_all; group(i).ctgs = [7]; group(i).legend = '100% Irr SchB';
        end
        
    end
    
    % Mode 6
    if ~isempty(strfind(fname_suffix,'bndAvsB'))
        i=0; clear group
        if isempty(gr_submode); gr_submode = 1; end
        
        if any(gr_submode == 1)
            i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = '60% SchA';
            i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = '60% SchB';
        end
        if any(gr_submode == 2)
            i=i+1; group(i) = gr_all; group(i).ctgs = [3]; group(i).legend = '50% SchA';
            i=i+1; group(i) = gr_all; group(i).ctgs = [4]; group(i).legend = '50% SchB';
        end
        if any(gr_submode == 3)
            i=i+1; group(i) = gr_all; group(i).ctgs = [5]; group(i).legend = 'TCent SchA';
            i=i+1; group(i) = gr_all; group(i).ctgs = [6]; group(i).legend = 'TCent SchB';
        end
    end
    
    % Mode 7 - highly redundant with mode 9 - delete?
    if ~isempty(strfind(fname_suffix,'bndctg60'))
        i=0; clear group
        if isempty(gr_submode); gr_submode = 1; end
        
        if any(gr_submode == 1)
            i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = '60% Ctg1';
            i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = '60% Ctg2';
        end
        
        if any(gr_submode == 2)
            i=i+1; group(i) = gr_all; group(i).ctgs = [3]; group(i).legend = '60% Ctg3';
            i=i+1; group(i) = gr_all; group(i).ctgs = [4]; group(i).legend = '60% Ctg4';
        end
        
        if any(gr_submode == 3)
            i=i+1; group(i) = gr_all; group(i).ctgs = [5]; group(i).legend = '60% Ctg5';
            i=i+1; group(i) = gr_all; group(i).ctgs = [6]; group(i).legend = '60% Ctg6';
        end
        
        if any(gr_submode == 4)
            i=i+1; group(i) = gr_all; group(i).ctgs = [7]; group(i).legend = '60% Ctg7';
            i=i+1; group(i) = gr_all; group(i).ctgs = [8]; group(i).legend = '60% Ctg8';
        end
        
        if any(gr_submode == 5)
            i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = '100% Ctg1';
            i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = '100% Ctg2';
        end
        
        if any(gr_submode == 6)
            i=i+1; group(i) = gr_all; group(i).ctgs = [3]; group(i).legend = '100% Ctg3';
            i=i+1; group(i) = gr_all; group(i).ctgs = [4]; group(i).legend = '100% Ctg4';
        end
        
        if any(gr_submode == 7)
            i=i+1; group(i) = gr_all; group(i).ctgs = [5]; group(i).legend = '100% Ctg5';
            i=i+1; group(i) = gr_all; group(i).ctgs = [6]; group(i).legend = '100% Ctg6';
        end
        
        if any(gr_submode == 8)
            i=i+1; group(i) = gr_all; group(i).ctgs = [7]; group(i).legend = '100% Ctg7';
            i=i+1; group(i) = gr_all; group(i).ctgs = [8]; group(i).legend = '100% Ctg8';
        end
    end
    
    % Mode 8
    if ~isempty(strfind(fname_suffix,'nonmtch_incgt'))
        i=0; clear group
        if isempty(gr_submode); gr_submode = 1; end
        
        if any(gr_submode == 1)
            i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = 'Non-match congruent';
            i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = 'Non-match incongruent';
        end
        
    end
    
    
    
    
    % Mode 9 - highly redundant with mode 9 - delete?
    if ~isempty(strfind(fname_suffix,'bndctgfull'))
        i=0; clear group
        if isempty(gr_submode); gr_submode = 1; end
        
        j=0;
        % 60%
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg1';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg2';
        end
        
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg3';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg4';
        end
        
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg5';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg6';
        end
        
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg7';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg8';
        end
        
        % 60%
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg1';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg2';
        end
        
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg3';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg4';
        end
        
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg5';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg6';
        end
        
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg7';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg8';
        end
        
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg1';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg2';
        end
        
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg3';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg4';
        end
        
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg5';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg6';
        end
        
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg7';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg8';
        end
        
        j=j+1;
        if any(gr_submode == j)
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '50% SchA';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '50% SchB';
        end
        
    end
    
    
    
    
end