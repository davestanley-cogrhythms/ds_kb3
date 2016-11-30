

function file_unit_names = fname_uname_merge (files_examined,units_examined)
%     file_unit_names = fname_uname_merge (files_examined,units_examined)

    file_unit_names = cellfunu( @(f,u2) cellfunu(@(u) [f ' ' u],u2), files_examined,units_examined); % Combine file names and unit names... God, how ugly...
end