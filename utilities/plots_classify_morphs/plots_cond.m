



% Analyzes image morphs by attempting to cluster using PCA/ICA
% % This is the method I end up using to classify test images!

%% Setup variables and load data.
    % This also activates automatic test plottin code

addpath(genpath(fullfile('..','funcs_classify_morphs')));
path_roy = getpath('path_roy');
path_conds = fullfile(path_roy,'morphs','trial_conds');

[conds, fnames, textdata] = load_conds(path_conds, 1);

