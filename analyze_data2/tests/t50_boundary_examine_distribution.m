



%% Testing plot for mode 5
i=5; figure; subplot(121); hist(morphA(ctgsetli(:,i)),[0:5:100]); xlim([-5 105]); title('SchA Morph Axis'); subplot(122); hist(morphB(ctgsetli(:,i)),[0:5:100]);  xlim([-5 105]); title('SchB Morph Axis'); 


%% Testing plot for mode 4
i = 8;
ind = schA; figure; subplot(121); hist(morphA(ctgsetli(:,i) & ind),[0:5:100]); xlim([-5 105]); title('SchA Morph Axis'); subplot(122); hist(morphB(ctgsetli(:,i) & ind),[0:5:100]);  xlim([-5 105]); title('Sch B Morph Axis'); legend('Sch A Trials');
ind = schB; figure; subplot(121); hist(morphA(ctgsetli(:,i) & ind),[0:5:100]); xlim([-5 105]); title('SchA Morph Axis'); subplot(122); hist(morphB(ctgsetli(:,i) & ind),[0:5:100]);  xlim([-5 105]); title('Sch B Morph Axis'); legend('Sch B Trials');

