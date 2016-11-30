



% runSGE.m is a script m-file to open matlabpool and run user app
% n and m are set at runtime as input to single_node_batch script
matlabpool open local        % use local configuration

% matmulExample(n, m)          % replace with your script or function m-file

currdir = pwd;
cd ..
setup_paths
setup_paths
cd (currdir);

identify_clipping(curr_stage);



