
function [gr, mask] = group_setup_FR2D_old
    % Group coordinates are a 4D matrix ordered based on relevance of
    % changes in firing rate significances


    
    % Cell selection preset
    ps_stage_any = [2 2 2 2];

    % Trial selection preset
    all_correct_trials = [1 2 3 4];

    % % Templates. Format is: AdA BdB AdB BdA for current and alternate stage, respectively. Format XdY means cue scheme X with significant difference across scheme Y.
    % % One decider
    % one_decider.criteria = [1 0 ps_nosig ps_stage_any]; one_decider.ctgs = [1 2];
    % one_decider.criteria = [0 1 ps_nosig ps_stage_any]; one_decider.ctgs = [3 4];
    % 
    % % Both deciders
    % both_deciders.criteria = [1 1 ps_nosig ps_stage_any]; both_deciders.ctgs = [1 2 3 4];

    % All possibilities Num right vs num wrong grouping criteria
    clear gr.all
    Ngrall = 4;
    sz = [sqrt(Ngrall) sqrt(Ngrall)];
    for i = 1:Ngrall
        [gr.all(i).coords(1), gr.all(i).coords(2)] = ind2sub(sz,i);
        gr.all(i).criteria = [coords2criteria(gr.all(i).coords) ps_stage_any];
        gr.all(i).ctgs = [all_correct_trials];
        [gr.all(i).right, gr.all(i).wrong] = coords2score(gr.all(i).coords);
        gr.all(i).diff = gr.all(i).right + gr.all(i).wrong;
        gr.all(i).netprefA = coords2preference(gr.all(i).coords);
        gr.all(i).metadata.sz = sz;
        gr.all(i).legend = '';
    end
    
    mask = [0 1; ...
            1 0];
    mask = logical(mask);

    % Scored by Right minus wRong
    %gr.all = normalize_diffs(gr.all,'diff');
    gr.rmr = score2grouping(gr.all,[gr.all.right]);
    % gr_netpref = score2grouping(gr.all,[gr.all.netprefA]); gr_netpref.ctgs = []
    
    
    % % All cells
    all_cells.criteria = [2 2 2 2 2 2 2 2];     all_cells.ctgs = [1 2 3 4];     all_cells.legend = 'All Cells';
    all_cells.coords = [1,1];
    gr.all_cells = all_cells;
    
    % % % % % % % % % % % % % % % % % % 
  
    
    % % % % % % % % % % % % % % % % % % 
    % % Original plots_preferred prototypes
    % % % % % % % % % % % % % % % % % % 
    % % A or B or Both
    clear one_or_both
    one_or_both.criteria(1,:) = [1 0 2 2 [2 2 2 2]]; one_or_both.ctgs(1,:) = [1 2 3 4]; one_or_both.legend = ('Preference');
    one_or_both.criteria(2,:) = [0 1 2 2 [2 2 2 2]]; one_or_both.ctgs(2,:) = [1 2 3 4];
    one_or_both.criteria(3,:) = [1 1 2 2 [2 2 2 2]]; one_or_both.ctgs(3,:) = [1 2 3 4];
    % 
    % % None all
    gr.none_all.criteria = [0 0 0 0 0 0 0 0]; gr.none_all.ctgs = [1 2 3 4]; gr.none_all.legend = ('Nothing in Sample nor Delay');
    gr.none_all.coords = [3 3];
    % 
    % % None curr
    gr.none_curr.criteria = [0 0 0 0 [2 2 2 2]]; gr.none_curr.ctgs = [1 2 3 4];
    gr.none_curr.coords = [3 3];
    % 
    

end



function crit = coords2criteria(coords)    

    use_jefferson_format = 1;   % Places Arel and Brel first, followed by A and B irrel
    % crit = [schemeA schemeB]
    % in teh graph, schemeA is columns, schemeB is rows
    crit = [coords2crit_scheme(coords(1)) coords2crit_scheme(coords(2))];
    
    
    % Rearrange to use Jefferson's format
    if use_jefferson_format
        crit = [crit(1) crit(3) crit(4) crit(2)];
            % Recall Jefferson's format goes:
            % SchACtg1/2    SchBCtg1/2  SchA_Nr_Ctg1/2  SchB_Nr_Ctg1/2
            %   AdA             BdB          BdA             AdB
            %  (It's convoluted as heck!)
    end
end


function crit = coords2crit_scheme(coords)
    switch coords
        case 1; crit = [1 2]; 
        case 2; crit = [0 2];
        otherwise
            fprintf('coords must be between 1 and 2. Exiting \n');
            crit = -1;
            return
    end
        
end


function [right, wrong] = coords2score(coords)

    right = sum(coords == 1);
    wrong = - sum(coords == 4); % Make wrongs a negative number for future
    
end


function [overall_prefA] = coords2preference(coords)

    prefA = sum(coords(:,1) == 1) + sum(coords(:,2) == 4);
    prefB = sum(coords(:,1) == 4) + sum(coords(:,2) == 1);
    
    overall_prefA = prefA - prefB;
    
end
