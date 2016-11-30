

function s = fix_freqband_format (s)

    % If only have 1 entry for freqband or timeband data, increase to 2
    % entries (so as to express it as a range)
    
    s.freqband_stats = conditional_clone(s.freqband_stats);
    s.freqband_perm = conditional_clone(s.freqband_perm);
    
    s.timeband_stats = conditional_clone(s.timeband_stats);
    s.timeband_perm = conditional_clone(s.timeband_perm);
    

end

function x= conditional_clone(x)
    % Convert from x=[var] to x=[var var]
    if length(x) == 1;
        x = repmat(x,[1,2]);
    end
    
end