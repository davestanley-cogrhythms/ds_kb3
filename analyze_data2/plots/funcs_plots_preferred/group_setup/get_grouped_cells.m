function [group_cells, num_cells] = get_grouped_cells(gr,scheme_preferred_all)


%     Criteria contains the sets of criteria we want in our cell subgroups.
%     Criteria are along columns, cell subgroups are added by adding new rows.
%     The resulting cell groups are specified in group_cells, wiht a new
%     column for each cell subgroup.

%     spc and spa are current and alternate preferred schemes (i.e. sample
%     and delay) for each neuron. The format of spc
%     is [Arel Airrel Brel Airrel], and likewise for spa for all cells (in rows).

    criteria_all = [gr.criteria gr.criteria_alt gr.criteria_sfc];
    %criteria_all = [gr.criteria gr.criteria_alt];

    %scheme_preferred_all = [spc spa sfcp];
    
    Ncells = size(scheme_preferred_all,1);
    Nsubgroups = size(criteria_all,1);
    Ncriteria = size(criteria_all,2);
    Ncriteria2 = size(scheme_preferred_all,2);
    
    if Ncriteria ~= Ncriteria2
        error('Error, mismatch in size of spc/spa and group.criteria.');
    end
    
    group_cells = true(Ncells,Nsubgroups);
    
    
    for i = 1:Nsubgroups
        for j = 1:Ncriteria
            curr_group = criteria_all(i,j);
            if curr_group < 2       % curr_group == 2 means we don't care whether the neuron was selective for this situation
                group_cells(:,i) = group_cells(:,i) & (scheme_preferred_all(:,j) == curr_group);
            end
        end
    end
    
    
    
    num_cells = sum(group_cells);

end

