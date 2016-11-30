
function myunits = lfppairs2unitpairs(mypairs,units2elect)

        % Build elec2units mapping (inverse mapping of units2elect)
        Nelects = max(mypairs(:));
        elec2units = get_elec2units(units2elect,Nelects);
        
        % Identify units corresponding to each of mypairs2
        myunits = [elec2units(mypairs(:,1)) elec2units(mypairs(:,2))];
        
end