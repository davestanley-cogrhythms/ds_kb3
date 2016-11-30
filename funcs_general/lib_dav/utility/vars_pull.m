function vars_pull(s,chosen_fields)

    if nargin < 2
        chosen_fields = fieldnames(s);
    end
    
    chosen_fields = chosen_fields(:)';

    for n = chosen_fields
        name = n{1};
        value = s.(name);
        assignin('caller',name,value);
    end
    
end