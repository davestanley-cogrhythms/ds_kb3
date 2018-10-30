


restoredefaultpath;

try
  rng('shuffle'); % - exception: rand in parfor does not use this new seed
  % rng(now);
catch
  %if any(regexp(version,'(19\d\d)|(200\d)|(2010)')) % if version before 2011a
  warning('on','MATLAB:RandStream:ActivatingLegacyGenerators')
  warning('on','MATLAB:RandStream:ReadingInactiveLegacyGeneratorState') 
  rand('seed',sum(100*clock));
end

% Runs the script containing the list of all the paths used by this
% software. This script returns variables identifying relevant paths

currpath = pwd;
addpath(currpath);
addpath(genpath(fullfile(currpath,'funcs_general')));

my_paths

% External paths
addpath(genpath(chronuxpath));              % Chronux
addpath(genpath(sfcpath));                  % Kyle's SFC code
addpath(genpath(path_src,'lib_dav'));       % Dave's library
addpath(genpath(path_src,'SigProc-Plott')); % SigProc-Plott code (make sure you're in dev branch!)


