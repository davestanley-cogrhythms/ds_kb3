
function classify_morphs

    % For each file, this maps the condition codes onto image ctgs for the test
    % images

    % Setup general stuff
    format compact
    format long g
    run_setdefaultfig
    addpath(genpath(fullfile('.','funcs_classify_morphs')));

    % Setup paths
    path_roy = getpath('path_roy');
    outpath = getpath('path_conditions');

    % First, classify the test images (under morphs/trial_conds/images)
    % as cat/dog/fat/thin using an SVM.
    % Returns the ctg of each test image.
    [schActg, schBctg, item_nums,files_all] = calc_morphs_knorm;

    % Sort the test image files based on item number
    [~,I] = sort(item_nums);
    schActg = schActg(I);
    schBctg = schBctg(I);
    item_nums = item_nums(I);
    files_all = files_all(I);

    % Load the conditions files (*.con, under morphs/trial_conds)
    path_conds = fullfile(path_roy,'morphs','trial_conds');
    [conds, fnames, textdata] = load_conds(path_conds);

    % The item numbers in each condition file 
    % correspond to the test images classified above. Thefefore, we can figure
    % out the test image type for each condition code.
    for i = 1:length(fnames)
        col=5; conds_to_test1{i} = map_conds2ctg(conds{i},schActg,schBctg,item_nums,col);

    end
    for i = 1:length(fnames)
        col=6; conds_to_test2{i} = map_conds2ctg(conds{i},schActg,schBctg,item_nums,col);
    end
    
    info_str = 'Computed by classify_morphs.m. conds_to_test1 and conds_to_test2 are for test1 and test2 images, respectively. Columns are as follows: Col 1 is condition code; Col 2 is 1 for cat, 0 for dog; Col 3 is 1 for thin, 0 for fat.';
    mkdir_dave(outpath);
    save(fullfile(outpath,'conditions.mat'),'fnames','conds_to_test1','conds_to_test2','textdata','info_str');
    
end


function conds2test_ctg = map_conds2ctg(conds,schActg,schBctg,item_nums,col)
    
    conds2test_ctg = zeros(size(conds,1),3);
    
    for i = 1:size(conds,1)
        
        ind = find(item_nums == conds(i,col));
        if ~isempty(ind); 
            conds2test_ctg(i,2:3) = [schActg(ind) schBctg(ind)];
        else
            if conds(i,col) ~= -4; warning('Item number in condition file not found'); end
            conds2test_ctg(i,2:3) = [NaN NaN];
        end
        
    end
    
    conds2test_ctg(:,1) = conds(:,1);
    
end


