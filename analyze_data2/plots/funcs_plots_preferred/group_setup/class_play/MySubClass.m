
classdef MySubClass < MySuperClass

    properties
        y
    end
    
    methods
        function B = MySubClass
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