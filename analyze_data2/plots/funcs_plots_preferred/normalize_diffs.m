function gr_all_out = normalize_diffs(gr_all,field)
    if nargin < 2
        field = 'diff';
    end

    gr_all_out = gr_all;
    for i = 1:length(gr_all)
        gr_all_out(i).(field) = 2*heaviside(gr_all(i).(field))-1;
    end

end