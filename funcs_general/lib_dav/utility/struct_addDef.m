function s = struct_addDef(s,fieldname,default_value)
    % If the field "fieldname" doesn't already exist in structure s, adds
    % it and assumes the default value. If defualt_value isn't specified,
    % adds an empty vector
    if nargin < 3
        default_value = [];
    end
    
    if ~isfield(s,fieldname); s.(fieldname) = default_value; end
end