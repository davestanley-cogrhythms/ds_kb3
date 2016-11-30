
%% 

%save('test_ntrials','currlfp_all','params');
load test_ntrials

N0 = [5:15];
M=100;

currlfp1 = currlfp_all(:,:,1);
currlfp2 = currlfp_all(:,:,2);

Cave_tr = [];
for i = 1:length(N0)
    Ntrials = N0(i);
    ind = round(unifrnd(1,size(currlfp1,2),Ntrials,M));
    
    Cave_temp = [];
    for j = 1:M
        [Cave,phi,S12ave,S1ave,S2ave,f]=coherencyc(currlfp1(:,ind(:,j)),currlfp2(:,ind(:,j)),params);
        Cave_temp = [Cave_temp, Cave(:)];
    end
    Cave_tr = [Cave_tr, mean(Cave_temp,2)];
end

figure; plot(Cave_tr)


%% 

