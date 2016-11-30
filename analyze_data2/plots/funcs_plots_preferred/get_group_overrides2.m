

function group = get_group_overrides2(ctgsetli_mode, N_criteria,xlims_desired0,ctgsetli_mode2)
    % Simplified version of original get_group_overrides
        % This just generates a group for each of the ctgs, and populates
        % the legend accordingly.
        % Doesn't fix the ordering (i.e. revers 50% and 100%) like previous
        % version did

    
    if ~exist('N_criteria','var'); N_criteria = 3; end      % Pick a random number
    if ~exist('ctgsetli_mode2','var'); ctgsetli_mode2 = 0; end      % Set to zero by defualt (ctgsetli modes 1-10)
    if exist('xlims_desired0','var');
        gr_all.xlims_desired = xlims_desired0;
    end
    
    group = []; % Default output

    % Set up group defaults
    gr_all = Grp;
    gr_all.criteria = [2*ones(1,N_criteria)];
    gr_all.criteria_alt = [];
    gr_all.criteria_sfc = [];
    gr_all.ctgs = 1;
    

    % Mode 0 - Default Ctgs - i.e. Jefferson
    if ctgsetli_mode == 0
        switch ctgsetli_mode2
            case {0,1}
                i=0; clear group
        %         i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Ctg1';
        %         i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Ctg2';
        %         i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Ctg3';
        %         i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Ctg4';

                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Cat';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Dog';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Goc';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Tad';

                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Cat irr';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Dog irr';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Goc irr';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Tad irr';

                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'SchA';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'SchB';

                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Ctg11';
            case 4
                i=0; clear group
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Cat';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Cat irr';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Dog';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Dog irr';

                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Goc';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Goc irr';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Tad ';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Tad irr';

                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'SchA';
                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'SchB';

                i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Ctg11';
        end
    end
    
    % Mode 1 - Switch trials
    if ctgsetli_mode == 1
        i=0; clear group
        i=i+1; group(i) = gr_all; group(i).ctgs = [5]; group(i).legend = 'Non-switch Trials';
        i=i+1; group(i) = gr_all; group(i).ctgs = [6]; group(i).legend = 'Switch Trials';
    end

    % Mode 2 - Reaction time
    if ctgsetli_mode == 2
        clear group
        i=0;
        i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = 'Fast RT';
        i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = 'Slow RT';            

    end

    % Mode 3 - reaction time2
    if ctgsetli_mode == 3
        clear group
        i=0;
        
        i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = 'Congruent Fast RT';
        i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = 'Congruent Slow RT';
        i=i+1; group(i) = gr_all; group(i).ctgs = [3]; group(i).legend = 'Incongruent Fast RT';
        i=i+1; group(i) = gr_all; group(i).ctgs = [4]; group(i).legend = 'Incongruent Slow RT';
           
    end



    % Mode 4 - boundary trials
    if ctgsetli_mode == 4
        i=0; clear group
        i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = '100% rel.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = '60% rel.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [3]; group(i).legend = '100% irr.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [4]; group(i).legend = '60% irr.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [5]; group(i).legend = '100% rel.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [6]; group(i).legend = '50% rel.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [7]; group(i).legend = '100% irr.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [8]; group(i).legend = '50% irr.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [9]; group(i).legend = '100% rel.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [10]; group(i).legend = '50 or 60% rel.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [11]; group(i).legend = '100% irr.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [12]; group(i).legend = '50 or 60% irr.';
    end
    
    % Mode 5 - Boundary Sch A Sch B
    if ctgsetli_mode == 5
        i=0; clear group
        i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = '100% C/D';
        i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = '50% C/D';
        i=i+1; group(i) = gr_all; group(i).ctgs = [3]; group(i).legend = '100% G/D';
        i=i+1; group(i) = gr_all; group(i).ctgs = [4]; group(i).legend = '50% G/D';
        i=i+1; group(i) = gr_all; group(i).ctgs = [5]; group(i).legend = '100% G/D irr.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [6]; group(i).legend = '50% G/D irr.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [7]; group(i).legend = '100% C/D irr.';
        i=i+1; group(i) = gr_all; group(i).ctgs = [8]; group(i).legend = '50% C/D irr.';
    end
    
    % Mode 6 - Sch A vs Sch B along boundaries
    if ctgsetli_mode == 6
        if ctgsetli_mode2 == 0
        i=0; clear group
        i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = '60% C/D';
        i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = '60% G/D';
        i=i+1; group(i) = gr_all; group(i).ctgs = [3]; group(i).legend = '50% C/D';
        i=i+1; group(i) = gr_all; group(i).ctgs = [4]; group(i).legend = '50% G/D';
        i=i+1; group(i) = gr_all; group(i).ctgs = [5]; group(i).legend = 'TCent C/D';
        i=i+1; group(i) = gr_all; group(i).ctgs = [6]; group(i).legend = 'TCent G/D';
        elseif any(ctgsetli_mode2 == [2,3])
            
            i=0; clear group
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Bnd2 SchA';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Bnd1 SchA';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Bnd2 SchB';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Bnd1 SchB';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Bnd2 SchA irr';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Bnd1 SchA irr';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Bnd2 SchB irr';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Bnd1 SchB irr';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Bnd1 SchA';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Outer A';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Bnd1 SchB';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Outer B';
            
        
        end
        
    end
    
    % Mode 7 - Boundry trials for Ctg1 vs Ctg2 - Highly redundant with mode 9 - delete?
    if ctgsetli_mode == 7
        if any(ctgsetli_mode2 == [0,1])
        i=0; clear group
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Cat';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Dog';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Goc';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Tad';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Cat irr';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Dog irr';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Goc irr';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Tad irr';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Cat';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Dog';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Goc';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Tad';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Cat irr';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Dog irr';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Goc irr';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Tad irr';
        elseif any(ctgsetli_mode2 == [2,3])
            i=0; clear group
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% C/D';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% C/D';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% G/D';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% G/D';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% G/D irr.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% G/D irr.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% C/D irr.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% C/D irr.';
        end
    end
    
    % Mode 8
    if ctgsetli_mode == 8
        if ctgsetli_mode2 == 0
        i=0; clear group
        i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = 'Non-match congruent';
        i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = 'Non-match incongruent';
        elseif any(ctgsetli_mode2 == [2,3])
            i=0; clear group
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '50% C/D';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% C/D';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '50% G/D';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% G/D';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '50% G/D irr.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% G/D irr.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '50% C/D irr.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% C/D irr.';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Cent SchA';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Outer SchA';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Cent SchB';
            i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = 'Outer SchB';
        end
    end
    
    
    
    
    % Mode 9 - Ctg1 vs Ctg2 for various boundary conditions - highly redundant with mode 7
    if ctgsetli_mode == 9
        i=0; clear group
        % 60%
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg1';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg2';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg3';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg4';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg5';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg6';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg7';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '60% Ctg8';
        % 80%
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg1';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg2';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg3';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg4';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg5';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg6';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg7';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '80% Ctg8';
        % 100%
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg1';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg2';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg3';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg4';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg5';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg6';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg7';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '100% Ctg8';
        % 100%
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '50% SchA';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '50% SchB';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '50% SchA & irr';
        i=i+1; group(i) = gr_all; group(i).ctgs = [i]; group(i).legend = '50% SchB & irr';
    end
    
    
end