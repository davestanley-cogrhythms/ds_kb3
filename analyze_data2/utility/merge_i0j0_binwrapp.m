
function merge_i0j0_binwrapp(varargin)
%     Inputs: sfc_mode,stage,files,i0,j0,embed_outname_i0j0
%     This script is compiled and then called by run_batch.sh
%     Alternatively, you can have run_batch.sh call run_analysis_wrapp.m
%          directly
% 
%     To compile, run:
%        restoredefaultpath
%        % BE SURE TO CHECK FOR CORE DUMPS!
%        cd ../../; setup_paths; cd analyze_data2/utility
%        %mcc -R -singleCompThread -mv -o run_analysis_binwrapp run_analysis_binwrapp.m -R -nodisplay   % (For serial)
%        mcc -mv -o merge_i0j0_binwrapp merge_i0j0_binwrapp.m -R -nodisplay			       % (For parallel)


    % matlabpool open local        % use local configuration

    % Don't run this, instad just make sure path is set up right when compiling
%     currdir = pwd;
%     cd ..
%     setup_paths
%     setup_paths
%     cd (currdir);
    
    % Convert all input arguments from strings to appropriate format.
    for i = 1:nargin
        if ischar(varargin{i})
            varargin{i} = eval(varargin{i});
        end
    end
    
    matlabpool local 8

    merge_i0j0(varargin{:})
    %run_analysis_simple(varargin{:})

    matlabpool close
    
end
