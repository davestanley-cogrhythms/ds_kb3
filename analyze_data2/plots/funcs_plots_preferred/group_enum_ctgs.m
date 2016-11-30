
function group = group_enum_ctgs(gr2,sz,fname_suffix,xlims_desired0)
    
    if nargin < 4
        xlims_desired0 = [0 80];
    end
    
    if length(sz) < 3; sz(3) = 1; end
    
    group(1:sz(3)) = gr2;
    for i = 1:length(group);
        group(i).ctgs = i;
        group(i).data_name = gr2(1).data_name;
        group(i).xlims_desired = xlims_desired0;
        
        if nargin > 2
            % Call get_group_overrides to get legend names
            N_criteria = 1; %Whatever..doesn't matter
            group_override = get_group_overrides(fname_suffix, N_criteria, 1:100);  % (Hack to return everything)
            if ~isempty(group_override)                     % It's one of our overrides - RT, switch, or boundary
                ind = [group.ctgs] == group(i).ctgs;        % Find appropriate ctg
                group(i).legend = group(ind).legend;
            else                                % It's the default scheme - A vs B
                group(i).legend = ['Ctg ' num2str(i)];
            end
        end
        
    end
    
    
    
end