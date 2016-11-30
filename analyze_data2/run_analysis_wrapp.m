
function run_analysis_wrapp(sfc_mode,stage,files,i0,j0,embed_outname_i0j0)
    % THIS FUNCTION IS DEPRECIATED - NOW USE setup_paths_n_run.m for
    % general script execution.
    % Called by run_batch.sh when not using compiled matlab.
    % Alternatively, you can compile and run run_analysis_binwrapp.m
    % matlabpool open local        % use local configuration

    currdir = pwd;
    cd ..
    setup_paths
    setup_paths
    cd (currdir);

    run_analysis(sfc_mode,stage,files,i0,j0,embed_outname_i0j0)

    % matlabpool close
    
end
