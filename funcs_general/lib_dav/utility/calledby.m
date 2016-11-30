
function out = calledby(depth)
% out = calledby(func_name)
% OVERVIEW
%     Queries the source of a function call, with variable depth, returning
%     the calling function's name.
% 
% FORMS
%     out = calledby
%     out = calledby(depth)
% 
% DESCRIPTION
%     out = calledby
%       Returns the name of the calling function
%     out = calledby(depth)
%       Returns the name of the calling function at depth depth
% 
% INPUTS
%     depth - Depth of calling function to search. 0=current function,
%     1=calling function (default), etc..
% 
% OUTPUTS
%     out - Name of calling function or queried function. Returns 'root' if
%     called from base workspace
% 
% EXAMPLES
%     Example 1
%     calledby                 % Returns error
%     calledby(0)              % Returns root
%     calledby(1)              % Returns error
%     
%     Example 2
%     % Drop this into a script called test1.m
%     function test1
%         calledby
%         fprintf('Should return root. Press any key... \n'); pause
%         
%         calledby(1)             
%         fprintf('Should return root. Press any key... \n'); pause
%         
%         calledby(0)             
%         fprintf('Should return test1. Press any key... \n'); pause
%         
%         test2
%     end
%     function test2
%         calledby
%         fprintf('Should return test1. Press any key... \n'); pause
%         
%         calledby(1)             
%         fprintf('Should return test1. Press any key... \n'); pause
%         
%         calledby(2)             
%         fprintf('Should return root. Press any key... \n'); pause
%         
%         calledby(0)             
%         fprintf('Should return test2. Press any key... \n'); pause
%         
%         %calledby(3)             
%         %fprintf('Should return error \n');
%     end
% 
%     SEE ALSO
%         is_calledby
%         http://www.mathworks.com/matlabcentral/fileexchange/51280-is-calledby-func-name-

    if nargin < 1; depth = [];end
    if isempty(depth); depth = 1; end;

    [ST,I] = dbstack;
    
    % Add dummy data for base workspace
    ST(end+1)=ST(end);
    ST(end).name = 'root';
    ST(end).file = 'root.m';
    ST(end).line = 1;
    
  
    if depth + 2 > length(ST)   % Add 2, one for callby and one to switch starting depth to 0
        error('Depth parameter exceeds actual depth of function calls.');
        %if depth + 2 > length(ST); depth = length(ST) - 2; end     % If depth is too deep, automatically rescale.
    end

    out = ST(depth+2).name;

    

end

