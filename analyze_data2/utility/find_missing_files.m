


%% Code to see missing files

% Suspected missing
cd /usr2/postdoc/stanleyd/Miller/roy/buffer/analyze_data7b/test2mode_52.700001_u0_permu/stage_2
d1 = dir('./*.mat');
n1 = {d1.name};

% Example of full file list
cd /usr2/postdoc/stanleyd/Miller/roy/buffer/analyze_data7b/mode_22.45171111_ee_bndABrel_mik_t1019_bs_permu_dbsPS_prtl/stage_2
d2 = dir('./*.mat');
n2 = {d2.name};

%% Search for missing files
clear missing_files 
missing_ind = [];
j=1;
for i = 1:length(n2)
    temp=strcmp(n1,n2{i});
    if sum(temp) == 0;
        missing_files{j} = n2{i};
        missing_ind = [missing_ind i];
        j=j+1;
    end
    
end







