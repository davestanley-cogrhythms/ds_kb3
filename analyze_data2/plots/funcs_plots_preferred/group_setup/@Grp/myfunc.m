

function group = myfunc (group, myfield, myfun)
% Perfomrs function myfun on the field myfield for all
% groups in the array.
% Example code:
% group2(1:3)=Grp;
% myf = @(x) x*2;
% group2 = group2.myfunc('xlims_desired',myf);
% This multiplies xlims_desired by 2x for all group2 in the array.

%     group.(myfield) = myfun(group.(myfield));
    for i = 1:length(group)
        group(i).(myfield) = myfun(group(i).(myfield));
    end
end