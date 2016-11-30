

function [all_cells, all_sens, all_sens4, single, single_sens] = group_setup_other

    % % % Other Groups
    
    % % All cells - Any Ctg
    gc = repmat(Grp,1,1);
    gc.criteria = [2 2 2 2 2 ];     gc.ctgs = [1 2 3 4];     gc.legend = ' ';
    gc.coords = [1,1];
    all_cells = gc;
    
    % % All cells - Ctg1v2
    gc = repmat(Grp,1,2);
    i = 1; gc(i).criteria = [2 2 2 2 2 ]; gc(i).ctgs = [1]; gc(i).coords = [1 1]; gc(i).legend = 'All Ctg1';
    i = 2; gc(i).criteria = [2 2 2 2 2 ]; gc(i).ctgs = [2]; gc(i).coords = [1 2]; gc(i).legend = 'All Ctg2';
    all_sens = gc;
    
    % % All cells - Ctg1v2v3v4
    gc = repmat(Grp,1,2);
    i = 1; gc(i).criteria = [2 2 2 2 2 ]; gc(i).ctgs = [1]; gc(i).coords = [1 1]; gc(i).legend = 'All Ctg1';
    i = 2; gc(i).criteria = [2 2 2 2 2 ]; gc(i).ctgs = [2]; gc(i).coords = [1 2]; gc(i).legend = 'All Ctg2';
    i = 3; gc(i).criteria = [2 2 2 2 2 ]; gc(i).ctgs = [3]; gc(i).coords = [1 2]; gc(i).legend = 'All Ctg3';
    i = 4; gc(i).criteria = [2 2 2 2 2 ]; gc(i).ctgs = [4]; gc(i).coords = [1 2]; gc(i).legend = 'All Ctg4';
    all_sens4 = gc;
    
    
    % % SchA comparison
    gc = repmat(Grp,1,2);
    i = 1; gc(i).criteria = [1 2 2 2 2 ]; gc(i).ctgs = [1 2]; gc(i).coords = [1 1]; gc(i).right = 1; gc(i).legend = 'Decider';
    i = 2; gc(i).criteria = [0 2 2 2 2 ]; gc(i).ctgs = [1 2]; gc(i).coords = [1 2]; gc(i).right = 0; gc(i).legend = 'Non-decider';
    single = gc;
    % 
    % % A sensitivity
    gc = repmat(Grp,1,2);
    i = 1; gc(i).criteria = [1 2 2 2 2 ]; gc(i).ctgs = [1]; gc(i).coords = [1 1]; gc(i).legend = 'Decider Ctg1';
    i = 2; gc(i).criteria = [1 2 2 2 2 ]; gc(i).ctgs = [2]; gc(i).coords = [1 2]; gc(i).legend = 'Decider Ctg2';
    single_sens = gc;
    


end