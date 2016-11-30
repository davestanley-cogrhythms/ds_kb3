

function [schActg, schBctg, item_nums,files_all] = calc_morphs_knorm
% [schActg, schBctg, item_nums,files_all] = calc_morphs_knorm
% Analyzes image morphs by attempting to cluster using PCA/ICA
% % This is the method I end up using to classify test images!
% Classification codes:
%     schActg: 1=cat, 0=dog
%     schBctg: 1=thin, 0=fat

    %% Setup
    
    plot_on = 0;
    plot_debug = 0;


    addpath(genpath(fullfile('..','funcs_classify_morphs')));
    path_roy = getpath('path_roy');
    path_images = fullfile(path_roy,'morphs','trial_conds','images');
    path_morphlines = fullfile(path_roy,'morphs','catdog psycho testing','morphlines');

    ctg1='rnd_ctg_1';
    ctg2='rnd_ctg_2';
    ctg3='rnd_ctg_3';
    ctg4='rnd_ctg_4';

    %% Load image groups

    [im1, files1] = load_images_groups(fullfile(path_images,ctg1));
    [im2, files2] = load_images_groups(fullfile(path_images,ctg2));
    [im3, files3] = load_images_groups(fullfile(path_images,ctg3));
    [im4, files4] = load_images_groups(fullfile(path_images,ctg4));

    imm = load_images_groups(fullfile(path_morphlines));

    Ntests = length(im1) + length(im2) + length(im3) + length(im4);
    ctg1_ind = build_logical(1:100,Ntests);
    ctg2_ind = build_logical(101:200,Ntests);
    ctg3_ind = build_logical(201:300,Ntests);
    ctg4_ind = build_logical(301:400,Ntests);
    schA_ind = ctg1_ind | ctg2_ind;
    schB_ind = ctg3_ind | ctg4_ind;

        % Save all images in one large cell array
    imalli = {im1{:} im2{:} im3{:} im4{:}};
    files_all = {files1{:},files2{:},files3{:},files4{:}};
    immi = imm;

    %% Reshape image groups

    im1 = reshape_image(im1);
    im2 = reshape_image(im2);
    im3 = reshape_image(im3);
    im4 = reshape_image(im4);
    imm = reshape_image(imm);


    %% Calculate feature

    im1 = calc_similarity(im1,imm);
    im2 = calc_similarity(im2,imm);
    im3 = calc_similarity(im3,imm);
    im4 = calc_similarity(im4,imm);

    %% Merge to images all

    imall = [im1 im2 im3 im4];



    %% PCA (works poorly)
    if plot_debug
        %imall = imall';     % We want "data" to be rows, 
        [coeff,sc,latent] = princomp(imall');    % 400 data points, 16384 dimensions; 

        score = sc;

        figure;
        hold on; plot(score(1:100,1),score(1:100,2),'b.')
        hold on; plot(score(101:200,1),score(101:200,2),'g.')
        hold on; plot(score(201:300,1),score(201:300,2),'r.')
        hold on; plot(score(301:400,1),score(301:400,2),'m.')
        end


    %% Plot morphline averages  (works well!)
        % Basicaly figures out which morphs in the morph line are associated
        % with ctg1. Then, it looks at the similarities of the
        % test image with each of these ctg1 morphs. Repeat for ctg2..ctg4
    morph_filenames = get_filelist(fullfile(path_morphlines),'*.bmp');
    ctg1_mInd = strfind(morph_filenames,'2t1'); ctg1_mInd = cellfun(@(x) ~isempty(x),ctg1_mInd);
    ctg2_mInd = strfind(morph_filenames,'5t4'); ctg2_mInd = cellfun(@(x) ~isempty(x),ctg2_mInd);
    ctg3_mInd = strfind(morph_filenames,'4t1'); ctg3_mInd = cellfun(@(x) ~isempty(x),ctg3_mInd);
    ctg4_mInd = strfind(morph_filenames,'5t2'); ctg4_mInd = cellfun(@(x) ~isempty(x),ctg4_mInd);

        % Takes the average of the similarities of each test image with
        % ctg1..ctg4. For each test image, takte the difference of its ctg1
        % similarity and ctg2 similarity; likewise for ctg3 and 4. This roughly
        % tells how catlike or doglike it is. Can do PCA if you want instead,
        % but just hard coding it seems to work better.
    score_custom = [ mean([imall(ctg1_mInd,:)' - imall(ctg2_mInd,:)'],2) mean([imall(ctg3_mInd,:)' - imall(ctg4_mInd,:)'],2) ];
    % temp = [ [imall(ctg1_mInd,:)' - imall(ctg2_mInd,:)'] [imall(ctg3_mInd,:)' - imall(ctg4_mInd,:)'] ]; [coeff,sc,latent] = princomp(temp); score_custom = sc(:,1:2);
    % temp = [ [imall(ctg1_mInd,:)'  imall(ctg2_mInd,:)'] [imall(ctg3_mInd,:)' imall(ctg4_mInd,:)'] ]; [coeff,sc,latent] = princomp(temp); score_custom = sc(:,1:2);

    
    if plot_on
        figure;
        hold on; plot(score_custom(1:100,1),score_custom(ctg1_ind,2),'b.')
        hold on; plot(score_custom(101:200,1),score_custom(ctg2_ind,2),'g.')
        hold on; plot(score_custom(201:300,1),score_custom(ctg3_ind,2),'r.')
        hold on; plot(score_custom(301:400,1),score_custom(ctg4_ind,2),'m.')
    end


    %% Train & plot SVM
    if plot_on
        figure; SVMSTructschA = svmtrain(score_custom([schA_ind],:),[ones(100,1);zeros(100,1)],'showplot','true');
        % schActg = svmclassify(SVMSTructschA,score_custom([schB_ind],:),'ShowPlot',true);

        figure; SVMSTructschB = svmtrain(score_custom([schB_ind],:),[ones(100,1);zeros(100,1)],'showplot','true');
        % schBctg = svmclassify(SVMSTructschB,score_custom([schA_ind],:),'ShowPlot',true);
    end
    %% Classify morphs with SVM
    SVMSTructschA = svmtrain(score_custom([schA_ind],:),[ones(100,1);zeros(100,1)],'showplot','false');
    schActg = svmclassify(SVMSTructschA,score_custom([1:end],:),'ShowPlot',false)';

    SVMSTructschB = svmtrain(score_custom([schB_ind],:),[ones(100,1);zeros(100,1)],'showplot','false');
    schBctg = svmclassify(SVMSTructschB,score_custom([1:end],:),'ShowPlot',false)';

    %% Fix classification of known images using knowledge from training
    schBctg(ctg3_ind) = 1;
    schBctg(ctg4_ind) = 0;

    %% Plot random images of classifications
    if plot_on
        % Choose category of classified images to plot
        chosen_ctg = 1;     % ctg1,2,3,4
        show_test = true;   % true for plot classified test images, false for training

        switch chosen_ctg
            case 1; inds = schActg == 1 & (show_test.*schB_ind | ~show_test.*schA_ind);
            case 2; inds = schActg == 0 & (show_test.*schB_ind | ~show_test.*schA_ind);
            case 3; inds = schBctg == 1 & (show_test.*schA_ind | ~show_test.*schB_ind);
            case 4; inds = schBctg == 0 & (show_test.*schA_ind | ~show_test.*schB_ind);
        end  

        inds = find(inds);

        N = 20;
        inds2 = inds(randperm(length(inds)));
        inds2 = inds2(1:N);

        figure;
        for i = 1:length(inds2)
            clf;
            imagesc(imalli{inds2(i)})
            pause
        end
    end

    %% Convert filenames to item numbers
    item_nums = filenames2item_numbers (files_all);

end