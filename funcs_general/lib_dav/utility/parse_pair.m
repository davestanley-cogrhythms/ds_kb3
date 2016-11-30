

function [val, myvarargin_remaining] = parse_pair(myvarargin,chosen_optional,default)
%     [val, myvarargin_remaining] = parse_pair(myvarargin,chosen_optional,default)
% 
%     Takes a set of input arguments, and searches for the name and value
%     pair specified by "chosen_optional". If this pair is found, the value
%     is returned in "val"; the cell array of input arguments is returned
%     minus this name value par in "myvarargin_remaining". Othewise, the 
%     default value is returned and myvarargin is returned unchanged in 
%     myvarargin_remaining.
%     
%     Example 1 - Something found
%     myvarargin = {1:10, sin(1:10), 'blah', 50, 'blah2', 80, cos(1:20)};
%     chosen_optional = 'blah';
%     default = -1;
%     [val, myvarargin_remaining] = parse_pair(myvarargin,chosen_optional,default)
%     
%     Example 2 - Something not found
%     myvarargin = {1:10, sin(1:10), 'blah', 50, 'blah2', 80, cos(1:20)};
%     chosen_optional = 'blah3';
%     default = -1;
%     [val, myvarargin_remaining] = parse_pair(myvarargin,chosen_optional,default)
    
    index = find(strcmp(myvarargin,chosen_optional));
    
    if isempty(index)
        val = default;
        myvarargin_remaining = myvarargin;
    else
        val = myvarargin{index+1};
        N = length(myvarargin);
        index2 = ~build_logical([index index+1],N);
        myvarargin_remaining = myvarargin(index2);
    end
    
    
end