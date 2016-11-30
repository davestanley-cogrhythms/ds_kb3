

% Submit this to the batch system to compile matlab, since the stupid
% deamon keeps killing my interactive attempts to compile. (ARGH!)
% Use the command
% qsub -l h_rt=2:00:00 single_node_batch "compile_matlab" localoutput


       restoredefaultpath
       restoredefaultpath
       restoredefaultpath
       currdir = pwd;
       cd ..;
       setup_paths;
       setup_paths;
       cd(currdir)
       addpath(genpath('funcs_supporting_local'));
       mcc -mv -o run_analysis_binwrapp run_analysis_binwrapp.m -R -nodisplay
