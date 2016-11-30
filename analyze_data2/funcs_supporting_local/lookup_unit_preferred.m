    
function pschemes = lookup_unit_preferred(filename,unitname,schemespath,schemesname)

    fileunit = [filename ' ' unitname];
    
    
    pref = load(fullfile(schemespath,schemesname));
    unitind = [];
    for i = 1:length(pref.file_unit_names)
        if strcmp(pref.file_unit_names{i},fileunit)
            unitind = [unitind i];
        end
    end
    
    if length(unitind) > 1
        fprintf('Uh oh - more than one matching unit found. Should debug at this point \n');
        return;
    elseif length(unitind) == 1
        pschemes = pref.sch_raw(unitind,:);
    else
        pschemes = [];
    end
    
    
    
end