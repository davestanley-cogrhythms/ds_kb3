
function group = group_enum_ctgs_diff(gr2,sz,fname_suffix,xlims_desired0)

    if nargin < 4
        xlims_desired0 = [0 180];
    end
    
    if length(sz) < 3; sz(3) = 1; end
    
    group(1:sz(3)) = gr2;
    for i = 1:length(group);
        group(i).ctgs = i;
        group(i).data_name = ['|\Delta ' gr2(1).data_name ' |'];
        group(i).xlims_desired = xlims_desired0;
        
        if nargin > 2
            % Call get_group_overrides to get legend names
            N_criteria = 1; %Whatever..doesn't matter
            gr_submode = i;
            group_override = get_group_overrides(fname_suffix, N_criteria, gr_submode);
            if ~isempty(group_override)         % It's one of our overrides - RT, switch, or boundary
                n1=group_override(1).legend;
                n2=group_override(2).legend;

                % Trim n1 Boundary
                n1 = strrep(n1,'morphs','');
                n1 = strrep(n1,'irr.','');
                n1 = strrep(n1,'rel.','');
                n1 = strrep(n1,'  ',''); % Remove double spaces
                n1 = strrep(n1,'   ','');
                
                % Trim n2 RT
                n2 = strrep(n2,'Incongruent','');
                n2 = strrep(n2,'Congruent','');
                group(i).legend = [n1,' vs ',n2];
            else                                % It's the default scheme - A vs B
                group(i).legend = ['Ctg ' num2str(i*2-1) ' vs ' num2str(i*2)];
            end
        end
        
    end
    
    
    
end