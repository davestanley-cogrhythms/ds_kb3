

% SETDEFAULTs1 sets the default structure values 
%    sout = SETDEFAULTS(s1, s2) merges s2 into
% s1, overwriting where necessary, and saves the
% result in sout

% Downloaded from
% http://stackoverflow.com/questions/38645/what-are-some-efficient-ways-to-combine-two-structures-in-matlab
% Credit: Jason S

function sout = struct_merge(s1,s2)
sout = s1;
for f = fieldnames(s2)'
    sout = setfield(sout,f{1},getfield(s2,f{1}));
end

% 
% function sout = setdefaults1(s,sdef)
% sout = sdef;
% for f = fieldnames(s)'
%     sout = setfield(sout,f{1},getfield(s,f{1}));
% end
