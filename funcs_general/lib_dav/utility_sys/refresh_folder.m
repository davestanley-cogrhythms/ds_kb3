


function refresh_folder (output_path)
    % Deletes folder contents then recreates

    fprintf(['Deleting all data in folder: \n' output_path '.\nPress any key to continue or CTRL-C to abort.\n']);
    pause

    if ispc
      system( ['del ' output_path '/*.*']);
      system( ['rmdir ' output_path]);
      system( ['mkdir ' output_path]);
    else
      system( ['rm ' output_path '/*']);
      system( ['rmdir ' output_path]);
      system( ['mkdir ' output_path]);
    end
    
    
end