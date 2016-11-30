function push_workspace()

c = getappdata(0, 'WORKSPACE_STACK');
if isempty(c)
    c = {};
end

% Grab workspace
w = evalin('caller', 'whos');
names = {w.name};
s = struct;
for i = 1:numel(w)
    s.(names{i}) = evalin('caller', names{i});
end

% Push it on the stack
c{end+1} = s;
setappdata(0, 'WORKSPACE_STACK', c);

