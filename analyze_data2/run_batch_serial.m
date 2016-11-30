
% Script for recreating actions of run_batch.sh, only for running in
% serial within matlab. Useful for filling in missing files.


%% Start parallel

% matlabpool open local 8       % use local configuration

%% Setup paths
% 
% currdir = pwd;
% cd ..
% setup_paths
% setup_paths
% cd (currdir);


%% File list and Ncells

rawpath = fullfile(getpath('path_lfp_sample_delay'));
fnm = dir(fullfile(rawpath,'*.mat'));
fnm = {fnm.name};
[~,fnmp] = cellfunu(@(f) fileparts(f),fnm);

% Ncells_arr = [6 8 4 2 5 1 6 5 5 4 7 1 3 5 6 9 9 9 8 9 6 6 6 4 7 4 6 4 5 3 3 4 6 7 4 6 1 5 5 5 13 15 15 11 8 13 7 9 13 7 11 9 11 9 4 4 11 8 9 12 10 11 11 5 4 3 6 5 6 9 7 11 11 7 9 9 8 7 2];
% Nlfp_pairs_arr=[3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 1 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3];
% Nlfp_pairs_arr=[21 28 21 28 28 6 28 28 28 21 28 28 28 15 15 45 45 66 28 28 28 28 28 28 21 28 28 15 15 6 28 28 28 28 15 28 1 28 28 28 91 120 120 120 91 120 91 91 120 28 120 120 120 120 120 120 120 120 120 120 120 120 120 28 15 120 78 120 120 91 91 120 120 120 91 120 120 120 28];
Nlfp_var=[7 8 7 8 8 4 8 8 8 7 8 8 8 6 6 10 10 12 8 8 8 8 8 8 7 8 8 6 6 4 8 8 8 8 6 8 2 8 8 8 14 16 16 16 14 16 14 14 16 8 16 16 16 16 16 16 16 16 16 16 16 16 16 8 6 16 13 16 16 14 14 16 16 16 14 16 16 16 8];


%% Setup inputs
stage = 2;

%%
if 0
    %% For all files
    for fi = 1:length(fnm)
        fi
        parfor j = 1:Ncells_arr(fi)
%         parfor j = 1:Nlfp_var(fi)
%         parfor j = 1:Nlfp_pairs_arr(fi)
            tic
            run_analysis(2.201410,stage,fi,[1,9,10,11],j,1);
%             run_analysis(41.63131,stage,fi,[1:4],j,1);
%             run_analysis(22.4413101,stage,fi,[1:12],j,1);
        
            toc
        end
    end
end

if 1
    %% For filling in missing files only
    % Can use this in conjunction with ./utility/find_missing_files
    % Useful when random batch jobs annoyingly fail.
    
    
    missing_ind = [4977  ]';
    
    parfor k = 1:length(missing_ind)
        k
        temp=0;
        missing_curr = missing_ind(k);
        for fi = 1:length(fnm)
            for j = 1:Nlfp_pairs_arr(fi)
                temp=temp+1;
                if temp == missing_curr
                    %keyboard
                    k
                    run_analysis(22.4718111,stage,fi,[1:16],j,1);
                end
            end
        end
    end
end
