
% Analyzes image morphs by attempting to cluster using PCA/ICA

%% Setup

addpath(genpath(fullfile('..','funcs_classify_morphs')));
path_roy = getpath('path_roy');
path_images = fullfile(path_roy,'morphs','trial_conds','images');
path_morphlines = fullfile(path_roy,'morphs','catdog psycho testing','morphlines');

ctg1='rnd_ctg_1';
ctg2='rnd_ctg_2';
ctg3='rnd_ctg_3';
ctg4='rnd_ctg_4';

%% Load image groups

im1 = load_images_groups(fullfile(path_images,ctg1));
im2 = load_images_groups(fullfile(path_images,ctg2));
im3 = load_images_groups(fullfile(path_images,ctg3));
im4 = load_images_groups(fullfile(path_images,ctg4));

%% Reshape image groups

im1 = reshape_image(im1);
im2 = reshape_image(im2);
im3 = reshape_image(im3);
im4 = reshape_image(im4);


%% Merge to images all

imall = [im1 im2 im3 im4];

%% Find pixels that are always zero
ind = all(imall == 0,2); % Indexes of pixesl that are all zero - we can dump these
imall = imall(~ind,:);


%% PCA
%imall = imall';     % We want "data" to be rows, 
[coeff,sc,latent] = princomp(imall');    % 400 data points, 16384 dimensions; 

score = sc;

figure;
hold on; plot(score(1:100,1),score(1:100,2),'b.')
hold on; plot(score(101:200,1),score(101:200,2),'g.')
hold on; plot(score(201:300,1),score(201:300,2),'r.')
hold on; plot(score(301:400,1),score(301:400,2),'m.')


% 
% 
% figure;
% hold on; plot(SCORE(1,1:100),SCORE(2,1:100),'b.')
% hold on; plot(SCORE(1,101:200),SCORE(2,101:200),'g.')
% hold on; plot(SCORE(1,201:300),SCORE(2,201:300),'r.')
% hold on; plot(SCORE(1,301:400),SCORE(2,301:400),'m.')



%% ICA (takes forever)
ICs=2;
tic
[z_ic A T mean_z] = myICA(imall,ICs);
toc


score = z_ic';



figure;
hold on; plot(score(1:100,1),score(1:100,2),'b.')
hold on; plot(score(101:200,1),score(101:200,2),'go')
hold on; plot(score(201:300,1),score(201:300,2),'rx')
hold on; plot(score(301:400,1),score(301:400,2),'m>')

