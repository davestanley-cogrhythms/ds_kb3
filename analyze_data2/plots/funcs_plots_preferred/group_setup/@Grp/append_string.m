

function gr = append_string(gr,myfield,mytext)
    %Adds text to an existing string in object
    
    gr.(myfield) = [gr.(myfield) mytext];
        
end

