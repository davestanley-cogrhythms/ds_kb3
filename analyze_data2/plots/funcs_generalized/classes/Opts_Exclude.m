
classdef Opts_Exclude
    % Class containing default settings
    properties
        % Unit exclusion for calculating bads
        exclude_clipping = 1;
        exclude_60 = 0;
        exclude_nans = 1;
        excludeL = 0; 
        excludeO = 0; 
        remove_dependent = 0;       % Remove dependent electrode pairs for FFC
    end
end
    

    

