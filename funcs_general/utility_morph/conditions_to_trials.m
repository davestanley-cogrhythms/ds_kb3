% Convert conditions to trials based on our newly determined morph
% percentages
function [trials] = conditions_to_trials(conditions)

    Nconditions = length(conditions);
    trials = zeros(8,length(conditions));

    for i = 1:Nconditions
        [sch, morphA, morphB] = condition2morph(conditions(i));
        
        [ctgnum, ctgnum_nr] = morph2ctg(sch, morphA, morphB);
        if ctgnum > 0; trials(ctgnum,i) = 1; end
        if ctgnum_nr > 0; trials(ctgnum_nr+4,i) = 1; end
    end
    
    trials = trials';
end


function [ctgnum, ctgnum_nr] = morph2ctg(sch, morphA, morphB)
    
if morphA == 50; A=-10;
elseif morphA <= 100 && morphA >= 60 A=0;
else A=1;
end

if morphB == 50; B=-10;
elseif morphB <= 100 && morphB >= 60 B=0;
else B=1;
end


if sch == 'A'
    ctgnum = A+1;
    ctgnum_nr = B+3;
elseif sch == 'B'
    ctgnum = B+3;
    ctgnum_nr = A+1;
else
    fprintf('Bad sch numb'); return
end

    
end