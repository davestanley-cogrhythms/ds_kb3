
function y = do_baseline_subtract(x,baseline_subtract)
    if baseline_subtract
        sz = size(x);
        y = x - repmat(mean(x,2),1,sz(2));
    else
        y = x;
    end
end
