
classdef MySuperClass

    properties
        x
    end
    
    methods
        function A = MySuperClass
        end
        

        function B = MySubClass(A)
            fprintf('test \n');
            B = MySubClass;
            B.x=1;
            B.y=2;
%             A = MySuperClass;
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

