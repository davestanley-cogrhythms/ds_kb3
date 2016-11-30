
function out = is_calledby(func_name)
% out = is_calledby(func_name)
% OVERVIEW
%     Tests identify of parent (calling) functions. Returns 1 if the active
%     function has a parent function calling it; else returns 0. Also tests
%     for specific parent names.
% 
% FORMS
%     [out] = is_calledby
%     [out] = is_calledby(func_name)
% 
% DESCRIPTION
%     [out] = is_calledby
%         Returns 1 if the active function has a parent; otherwise
%         returns 0.
%     [out] = is_calledby(func_name)
%         Returns 1 if the active function's parent matches func_name;
%         otherwise returns 0.
% 
% INPUTS
%     func_name - string indicating the name of a function. Set
%     func_name='*' to return 1 for being called by any function.
% 
% OUTPUTS
%     out - returns 1 or 0
% 
% EXAMPLES
%     Example 1
%     is_calledby                 % Returns 0
%     
%     Example 2
%     % Drop this into a script called test1.m
%     function test1
%         is_calledby             % Returns 0
%         fprintf('Should return 0 \n');
% 
%         is_calledby('*')        % Returns 0
%         fprintf('Should return 0 \n');
%         test2
%     end
%     function test2
%         is_calledby             % Returns 1
%         fprintf('Should return 1 \n');
% 
%         is_calledby('*')        % Returns 1
%         fprintf('Should return 1 \n');
%         
%         is_calledby('asdfasdf') % Returns 0
%         fprintf('Should return 0 \n');
%         
%         is_calledby('test1')    % Returns 1
%         fprintf('Should return 1 \n');
%         
%     end
% 
%     SEE ALSO
%         calledby
%         http://www.mathworks.com/matlabcentral/fileexchange/51281-calledby-depth-

if nargin < 1; func_name = [];end

[ST,I] = dbstack;

if ~(isempty(func_name) || strcmp(func_name,'*'))
    if length(ST) < 3
        out=0;
    else
        out = strcmp(ST(3).name,func_name);
    end
else
    out = length(ST) > 2;
end


end
