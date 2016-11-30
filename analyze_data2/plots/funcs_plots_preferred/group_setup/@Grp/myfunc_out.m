

function out = myfunc_out(group, myfunc)


    out_temp = myfunc(group(1));
    out = repmat(out_temp,size(out_temp));
    
    if length(group) > 1
        for i = 2:length(group)
            out(i) = myfunc(group(i));
        end
    end
    
end