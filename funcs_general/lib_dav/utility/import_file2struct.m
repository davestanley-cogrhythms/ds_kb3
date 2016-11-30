function s = import_file2struct(s,fname,fieldname, force_reload)

    if ~isfield(s,fieldname) || force_reload
        s.(fieldname) = load (fname);
    end

end