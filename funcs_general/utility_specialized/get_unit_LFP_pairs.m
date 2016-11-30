
function [unit_LFP_pairs unit_LFP_metadata] = get_unit_LFP_pairs(unit_names,lfp_names)
    % Returns a 3xN matrix
    % First row just numbers 1:N
    % Second row contains the LFP number of the corresponding unit in
    % "unit_names"
    % The third row (optional) gives the unit number as defined in unit_names:
    % For example, if PFC_Elect_5_unit_3 is the 6th unit listed in unit_names
    % then the corresponding column in unit_LFP_pairs would be [6; 5; 3].

    Nunits = length(unit_names);
    Nelectrodes = length(lfp_names);
    unit_LFP_pairs = zeros(2,Nunits);
    unit_LFP_pairs(1,:) = 1:Nunits; % Number 1:N
    
    for i = 1:Nunits
        LFP_index = -1;
        for j = 1:Nelectrodes
            %unit_names{i}
            %lfp_names{j}
            if ~isempty(strfind(unit_names{i},lfp_names{j}))
                LFP_index = j;
            end
        end
        
        if LFP_index == -1
            fprintf('Error: corresponding electrode not found for unit %s \n', unit_names{i});
            break;
        end

        unit_LFP_pairs(2,i) = LFP_index;
        
        
        % % Get unit number
        unit_index = strfind(unit_names{i},'_unit_') + length('_unit_');
        unit_num = str2num(unit_names{i}(unit_index:end));
        
        % % Dump everything into a metadata structure
        unit_LFP_metadata{i}.unit_str = unit_names{i};
        unit_LFP_metadata{i}.unit_num = unit_num;
        unit_LFP_metadata{i}.LFP_str = lfp_names{LFP_index};
        unit_LFP_metadata{i}.LFP_index = LFP_index;
    end
end
