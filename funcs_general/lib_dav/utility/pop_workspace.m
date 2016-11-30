
function pop_workspace(blank_slate)
    % If blank slate is true, clear the current workspace before importing
    % from the stack. If it's false, leave the current workspace as-is, and
    % import new data overtop.

if nargin < 1; blank_slate = true; end

% Pop last workspace off stack
c = getappdata(0, 'WORKSPACE_STACK');
if isempty(c)
    warning('Nothing on workspace stack');
    return;
end
s = c{end};
c(end) = [];
setappdata(0, 'WORKSPACE_STACK', c);

% Do this if you want a blank slate for your workspace
if blank_slate
    evalin('caller', 'clear');
end

% Stick vars back in caller's workspace
names = fieldnames(s);
for i = 1:numel(names)
    assignin('caller', names{i}, s.(names{i}));
end