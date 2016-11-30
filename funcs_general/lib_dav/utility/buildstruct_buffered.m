function sworkspace = buildstruct_buffered(sworkspace, buffer_path, fieldname, myfunc, force_recalc, force_reload)
% Evaluates myfunc and saves the output in the structure sworkspace as a variable named
% fieldname
% myfunc takes no arguments and returns a struct


    %  Force_reload
        %  0 - don't force load
        %  1 - force load
        % -1 - force to not reload, even if structure is missing (only use
        %      in special cases

    filename = [fullfile(buffer_path,fieldname) '.mat'];


    if ~exist(filename,'file') || force_recalc
        sout = myfunc();
        %eval([fieldname ' = s;']);
        mkdir_dave(buffer_path);
        save(filename,'sout','-v7.3');

    %     % Pull adat fields into current workspace and save these variables
    %     s_fieldnames = fieldnames(s);
    %     vars_pull(s);
    %     clear s
    %     mkdir_dave(buffer_path);
    %     save(filename,s_fieldnames{:});
    
        sworkspace.(fieldname) = sout;      % Now we reload automatically after recalculating.

    end

    % Load file and save variables in workspace, if not there already
    if (~isfield(sworkspace,fieldname) || force_reload == 1) && force_reload ~= -1
        load (filename);
        sworkspace.(fieldname) = sout;
    end


    %sworkspace = import_file2struct(sworkspace,filename,fieldname, reload);


end



