

function index = classify_scheme(sch,morphA,morphB,chosen_sch)

    exclude_boundary_cases = 1;
    
    if exclude_boundary_cases
        index = (sch == chosen_sch & round(morphA) ~= 50 & round(morphB) ~= 50);
    else
        index = (sch == chosen_sch);
    end

end