


%% Load MD
filelist = get_filelist(fullfile(getpath('path_metadata'),'sansunits'));
file_range = 1:length(filelist);
md_all = cellfun(@(filename) load (fullfile(getpath('path_metadata'),'sansunits',filename)),filelist(file_range));


%% Plot Ntrials by morph percentage
% Calculates the percentage of correct trials at each boundary level for
% all files

% For each file...
histAa = [];
histBa = [];
tic
for i = 1:length(md_all)
    i
    md = md_all(i);
    
    [sch, morphA, morphB] = arrayfunu(@condition2morph,md.each_trial_condition);
    schA = cellfun(@(x) x == 'A', sch);
    morphA = cell2mat(morphA);
    morphB = cell2mat(morphB);
    correct = md.each_trial_response == 0;
    sot = md.sample_on_trials(:);
    
    mA = morphA(schA == 1 & correct & sot);
    mB = morphB(schA == 0 & correct & sot);
    
    if sum(correct) ~= sum(sot & correct); warning('Correct trials that are not sample on trials?'); end
    
    morph_range=0:10:100;
    [histA,bins] = hist(mA,morph_range); if any(bins ~= morph_range); warning('Bin mismatch'); end
    [histB,bins] = hist(mB,morph_range); if any(bins ~= morph_range); warning('Bin mismatch'); end
    
    histAa = [histAa histA(:)];
    histBa = [histBa histB(:)];
end
toc


% Plot

figure; bar(morph_range,histAa); xlabel('Morph Percentage'); ylabel('Ntrials'); title('All files - Scheme A'); xlim([-10 110]);legend('File1','File2','...')
figure; bar(morph_range,histBa); xlabel('Morph Percentage'); ylabel('Ntrials'); title('All files - Scheme B'); xlim([-10 110]);legend('File1','File2','...')

%% Plot Ntrials by Ctg and Sch
% % Calculates the percentage of correct trials assocaited with Scheme A/B
% and Ctg 1,2,3,4

% For each file...
hScha = [];
hCtga = [];



tic
for i = 1:length(md_all)
    i
    md = md_all(i);
    
    [sch, morphA, morphB] = arrayfunu(@condition2morph,md.each_trial_condition);
    schA = cellfun(@(x) x == 'A', sch);
    morphA = cell2mat(morphA);
    morphB = cell2mat(morphB);
    
    correct = md.each_trial_response == 0;
    sot = md.sample_on_trials(:);
    bad_boundary = (morphA == 50 | morphB == 50);
    
    ctg1 = get_good_samples(md.ctg1_trials,md);
    ctg2 = get_good_samples(md.ctg2_trials,md);
    ctg3 = get_good_samples(md.ctg3_trials,md);
    ctg4 = get_good_samples(md.ctg4_trials,md);
    
    
    % Calculate Sch distribution
    hSch = hist(schA( correct & sot & ~bad_boundary),[0,1]); hSch=fliplr(hSch);
    hSch2 = [sum(ctg1 | ctg2) sum(ctg3 | ctg4)];
    if any(hSch ~= hSch2); warning('Error estimating Schemes'); end
    
    ind = correct & sot & ~bad_boundary;
    hCtg = [sum(schA(ind) & morphA(ind) > 50) ...
        sum(schA(ind) & morphA(ind) < 50) ...
        sum(~schA(ind) & morphB(ind) > 50)  ...
        sum(~schA(ind) & morphB(ind) < 50) ]; 
    hCtg2 = [sum(ctg1),sum(ctg2),sum(ctg3),sum(ctg4)];
    if any(hCtg ~= hCtg2); warning('Error estimating Ctgs'); end


    hScha = [hScha hSch(:)];
    hCtga = [hCtga hCtg(:)];
end
toc



figure; bar(hScha); xlabel('Scheme'); set(gca,'XTickLabel',{'A','B'}); ylabel('Ntrials'); title('All files - Scheme'); xlim([0 3]);legend('File1','File2','...')
figure; bar(hCtga); xlabel('Ctg'); ylabel('Ntrials'); title('All files - Scheme B'); title('All files - Ctg'); xlim([0 5]);legend('File1','File2','...')

