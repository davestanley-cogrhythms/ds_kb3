

function [unum] = fnum2unum1D(n,fnames,funames1D)

    chosen_fnames = fnames(n);
    
    out_arr = [];
    out_arr = zeros(length(funames1D),length(n));
    for i = 1:length(chosen_fnames)
        [temp1, currfname] = fileparts(chosen_fnames{i});
        out = cellfunu(@(x) strfind(x,currfname),funames1D);
        out = cellfun(@(x) ~isempty(x),out);
        out_arr(:,i) = out(:);
    end
    
    unum = any(out_arr,2);

end