

function group = arrayload (group, myfield, mydata)
% Loads mydata into myfield for all groups in array

    for i = 1:length(group)
        group(i).(myfield) = mydata;
    end
end