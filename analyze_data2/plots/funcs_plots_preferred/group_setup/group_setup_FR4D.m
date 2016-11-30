
function [gr, mask] = group_setup_FR4D(gr)
    % Group coordinates are a 4D matrix ordered based on relevance of
    % changes in firing rate significances

    N_ctgs_base = get_Nctgs_base;
    N_ctgs_extra = get_Nctgs_extra;

    
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
    Ngrall = 16;
    gr.all = repmat(Grp,1,Ngrall);
    sz = [sqrt(Ngrall) sqrt(Ngrall)];
    for i = 1:Ngrall
        [gr.all(i).coords(1), gr.all(i).coords(2)] = ind2sub(sz,i);
        gr.all(i).criteria = [coords2criteria(gr.all(i).coords, N_ctgs_extra)];
        gr.all(i).ctgs = [all_correct_trials];
        [gr.all(i).right, gr.all(i).wrong] = coords2score(gr.all(i).coords);
        gr.all(i).diff = gr.all(i).right + gr.all(i).wrong;
        gr.all(i).netprefA = coords2preference(gr.all(i).coords);
        gr.all(i).metadata.sz = sz;
        gr.all(i).legend = '';
    end
    
    mask =  [1 1 1 1; ...
            1 1 1 1; ...
            1 1 1 1; ...
            1 1 1 1; ];
    mask =  [0 0 1 0; ...
             0 0 1 0; ...
             1 1 1 1; ...
             0 0 1 0; ];
    mask =  [1 1 1 1; ...
            1 0 0 0; ...
            1 0 0 0; ...
            1 0 0 0; ];
%     mask =  [0 0 1 1; ...
%             0 0 1 1; ...
%             1 1 0 0; ...
%             1 1 0 0; ];
    mask = logical(mask);

end


function crit = coords2criteria(coords, N_ctgs_extra)    

    use_jefferson_format = 1;   % Places Arel and Brel first, followed by A and B irrel
    N_pairs_extra = round(N_ctgs_extra / 2);
    
    % crit = [schemeA schemeB]
    % in teh graph, schemeA is columns, schemeB is rows
    crit = [coords2crit_scheme(coords(1)) coords2crit_scheme(coords(2)) 2*ones(1,N_pairs_extra)];
    
    % Rearrange to use Jefferson's format
    if use_jefferson_format
        crit(1:4) = [crit(1) crit(3) crit(4) crit(2)];
            % Recall Jefferson's format goes:
            % SchACtg1/2    SchBCtg1/2  SchA_Nr_Ctg1/2  SchB_Nr_Ctg1/2
            %   AdA             BdB          BdA             AdB
            %  (It's convoluted as heck!)
    end
    
end

function crit = coords2crit_scheme(coords)
    switch coords
        case 1; crit = [1 0]; 
        case 2; crit = [1 1];
        case 3; crit = [0 0];
        case 4; crit = [0 1];
        otherwise
            fprintf('coords must be between 1 and 4. Exiting \n');
            crit = -1;
            return
    end
        
end


function [right, wrong] = coords2score(coords)

    right = 1*sum(coords == 1) + 0*sum(coords == 2);
    wrong = - sum(coords == 4); % Make wrongs a negative number for future
    
end


function [overall_prefA] = coords2preference(coords)

    prefA = sum(coords(:,1) == 1) + sum(coords(:,2) == 4);
    prefB = sum(coords(:,1) == 4) + sum(coords(:,2) == 1);
    
    overall_prefA = prefA - prefB;
    
end
