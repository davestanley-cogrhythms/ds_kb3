

function gr = group_enumerate_sens(gr,mask,symmetry_mode)
    % Enumerates original group structure to test across various
    % sensitivity cases
%     symmetry_mode = 1 - A
%     symmetry_mode = 2 - B
%     symmetry_mode = 3 - Group by preference
%     symmetry_mode = 4 - Test symmetry preferred
%     symmetry_mode = 5 - Test symmetry nonpreferred

    
    if nargin < 2
        symmetry_mode = 3;
    end

    % Scored by Scheme A vs Schmee B
    curr_score = 0;
    
    % Group testing sensitivity to image category
    i=1;gc(i) = group_tri(gr.all,mask,1);       gc(i).ctgs = repmat(1,size(gc(1).coords,1),1);  gc(i).legend = 'Preferred Ctg1';
    i=2;gc(i) = group_tri(gr.all,mask,1);       gc(i).ctgs = repmat(2,size(gc(1).coords,1),1);  gc(i).legend = 'Preferred Ctg2';
    gc = apply_symmetry(gc,gr.all,symmetry_mode);
    gr.sens_ctg = gc;
    
    % Group testing sensitivity to image category - non-relevant
    i=1;gc(i) = group_tri(gr.all,mask,1);       gc(i).ctgs = repmat(3,size(gc(1).coords,1),1);  gc(i).legend = 'Non-Preferred Ctg1';
    i=2;gc(i) = group_tri(gr.all,mask,1);       gc(i).ctgs = repmat(4,size(gc(1).coords,1),1);  gc(i).legend = 'Non-Preferred Ctg2';
    gc = apply_symmetry(gc,gr.all,symmetry_mode);
    gr.sens_ctg_np = gc;
    
    % Group testing sensitivity to image category - irrelevant
    i=1;gc(i) = group_tri(gr.all,mask,1);       gc(i).ctgs = repmat(7,size(gc(1).coords,1),1);  gc(i).legend = 'Preferred Irrel Ctg1';
    i=2;gc(i) = group_tri(gr.all,mask,1);       gc(i).ctgs = repmat(8,size(gc(1).coords,1),1);  gc(i).legend = 'Preferred Irrel Ctg2';
    gc = apply_symmetry(gc,gr.all,symmetry_mode);
    gr.sens_ctg_ir = gc;
    
    % Sensitivity testing sensitivity to categorization scheme
    i=1;gc(i) = group_tri(gr.all,mask,1);       gc(i).ctgs = repmat([1 2],size(gc(1).coords,1),1);  gc(i).legend = 'Preferred Sch';
    i=2;gc(i) = group_tri(gr.all,mask,1);       gc(i).ctgs = repmat([3 4],size(gc(1).coords,1),1);  gc(i).legend = 'Non-Preferred Sch';
    gc = apply_symmetry(gc,gr.all,symmetry_mode);
    gr.sens_sch = gc;
    
    % Sensitivity to irrelevant categorization scheme
    gc = gr.sens_sch;
    swap_map = [1 2 3 4; [5 6 7 8] ];
    swap_map_pref = unique(round(swap_map/2)','rows')'; % Maps the ctg indices onto pref indicies (i.e. src -> spc). For example 1-4 maps onto 1-2; etc.
    %gc = apply_symmetry(gc,gr.all,symmetry_mode);
    gc = group_swap.crit (gc,swap_map_pref);
    gc = group_swap.ctgs (gc,swap_map);
    gr.sens_sch_ir = gc;
    
    % Scheme A preferring cells in in Sch A vs Scheme B preferring cells in Sch B
    i=1;gc(i) = group_tri(gr.all,mask,1);       gc(i).ctgs = repmat([1 2],size(gc(1).coords,1),1);  gc(i).legend = 'Sch A Cells';
    i=2;gc(i) = group_tri(gr.all,mask,0);       gc(i).ctgs = repmat([3 4],size(gc(1).coords,1),1);  gc(i).legend = 'Sch B Cells';
    %gc = apply_symmetry(gc,gr.all,symmetry_mode);
    gr.AvB_cells= gc;
    
    % All decider cells in Sch A vs Sch B
    i=1;gc(i) = group_merge(gr.all(1),gr.all(2),gr.all(3)); gc(i).ctgs = repmat([1 2],size(gc(1).coords,1),1); gc(i).legend = 'Sch A Ctgs';
    i=2;gc(i) = group_merge(gr.all(1),gr.all(2),gr.all(3)); gc(i).ctgs = repmat([3 4],size(gc(1).coords,1),1); gc(i).legend = 'Sch B Ctgs';
    gr.AvB_ctgs = gc;
    
    % All decider cells in Sch A vs Sch B
    clear gc
    gc(1:2) = Grp;
    gc(1).ctgs = 9; gc(1).legend='Sch A';
    gc(2).ctgs = 10; gc(2).legend='Sch B';
    gr.AvB_everything=gc;
    
end




function gout = group_tri(group,mask,is_upper)
%     isupper = 1/0 for upper/lower triangle

    if nargin < 3
        is_upper = 1;
    end
    
    sz = group(1).metadata.sz;
    
    if is_upper tri = triu(true(sz),1);
    else tri = ~triu(true(sz),0);
    end
    tri=tri(:); mask=mask(:);
    ind = tri & mask;
    Ngout = sum(ind);
    gout = group_merge(group(ind));
    
    %gout.ctgs = repmat([1],Ngout,1);

    
end




function gout = apply_symmetry(group,gr_all, symmetry_mode)
%     symmetry_mode = 1 - A
%     symmetry_mode = 2 - B
%     symmetry_mode = 3 - Group by preference
%     symmetry_mode = 4 - Test symmetry preferred
%     symmetry_mode = 5 - Test symmetry unpreferred

    grA = group;
    grB = group_complement(grA,gr_all); % "B" complement of A

    switch symmetry_mode
        case 1; gout = grA;
        case 2; gout = grB;
        case 3
            for i = 1:length(grA)
                gout(i) = group_merge([grA(i), grB(i)]);
            end
        case 4; gout = [grA(1), grB(1)];
        case 5; gout = [grA(2), grB(2)];
    end
end
    