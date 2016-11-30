
function out = is_gpuArray(x)
    S = whos('x');
    out = strcmp(S.class,'gpuArray');
    
end