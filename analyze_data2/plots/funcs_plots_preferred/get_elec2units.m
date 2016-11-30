
function elec2units = get_elec2units(units2elect,Nelects)

        % List all available units
        units_numbers = [1:length(units2elect)]';
        
        % Build elec2units mapping (inverse mapping of units2elect)
        elec2units = Inf*ones(Nelects,1);  % Electrodes that dont have a unit are Infs
        elec2units(units2elect) = units_numbers;
        
end