
%%
% This is a script I wrote to search through my saved .mat files in the
% matfiles folder and identify mat files that are run using too low a value
% of Nperms (in this case, it was Nperms=2. ARGhhh!)

clearvars -except wrkspc_buffer
matfiles_path = getpath('path_matfiles_store');

n = dir(fullfile(matfiles_path,'per*'));
n={n.name};

uls = ones(1,length(n))*NaN;
uld = ones(1,length(n))*NaN;
for i = 1:length(n)
    load(fullfile(matfiles_path,n{i}));
    if isfield(sout,'stage2'); uls(i) = length(unique(sout.stage2.pvals_Cave(:)));end
    if isfield(sout,'stage3'); uld(i) = length(unique(sout.stage3.pvals_Cave(:)));end
end

figure;
plot(uls,'b');
hold on; plot(uld','r>'); legend('unique p vals for sample','unique p vals for delay');
xlabel('file #'); ylabel('N unique p vals');

% Work with uld, since it seems to be more reliable
ind = find(uld > 0 & uld < 80);
files_with_nperms_2 = n(ind)

% Missing stage2 (just for fun...)
ind2 = find(isnan(uls));
files_missing_stage2= n(ind2)

