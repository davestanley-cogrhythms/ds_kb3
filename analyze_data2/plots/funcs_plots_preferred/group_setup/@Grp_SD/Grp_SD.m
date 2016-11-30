
classdef Grp_SD < Grp
    % Group class that distinguishes between current and alternate
    % stage (i.e. sample and delay)

    properties
        
        % Minimal properties that must be specified (in instantiation)
        
        % Additional properties
        coords_alt
        right_alt
        wrong_alt
        diff_alt
        score_alt
        
        % Populated data
        
    end
    
    methods
        
        function gr = Grp_SD
            
        end
        
        function output = inheritObj(output,input)
            C = metaclass(input);
            P = C.Properties;
            for k = 1:length(P)
                if ~P{k}.Dependent
                    output.(P{k}.Name) = input.(P{k}.Name);
                end
            end
        end
    end

end


