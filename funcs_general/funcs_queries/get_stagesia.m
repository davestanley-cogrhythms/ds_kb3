


function [boundsia] = get_stagesia(stage,sample_on_ind)
    % Pull out boundaries of a given stage
    % Boundary outputs are in indices!
    % Boundary outputs are absolute relative to the file's lfp indices
    % Note that, by convention, all data in a given
    % stage should be defined as >= boundsir(1) and < boundsir(2).
    
    % Meaning of stage:
    % Notation is derived from get_stagesir. See get_stagesir.m
    
    boundsir = get_stagesir(stage);
    boundsia = repmat(boundsir(:)',length(sample_on_ind),1) + repmat(sample_on_ind(:),1,2);

end