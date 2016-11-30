

function [fnum, unum, curr_funame] = unum2fnum(n,fnames,funames)
    
    funames_1d = cat(2,funames{:});
    curr_funame = funames_1d{n};
    
    findex = zeros(1,length(fnames));
    for i = 1:length(fnames)
        out = strfind(curr_funame,fnames{i});
        findex(i) = ~isempty(out);
    end
    fnum = find(findex);
    if length(fnum) >1 fprintf('Error, duplicate filenames. \n'); return; end
    
    
    unum = find(strcmp(curr_funame,funames{fnum}));
    
    curr_funame = funames_1d{n};
    
    %Sanity check:
    test_funame = [funames{fnum}{unum}];
    
    if ~strcmp(curr_funame,test_funame)
        fprintf('Sanity check failed - something wrong! \n');
        return;
    end

end