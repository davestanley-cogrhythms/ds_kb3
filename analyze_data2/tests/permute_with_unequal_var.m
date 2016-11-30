


%% Standard permutation test

m=4;
n=100;
Nperms = 1500;
xstd=10;

x=xstd*randn(m,1);
y=randn(n,1);

% Permutation test
comb = [x;y];
xp=[];
yp=[];
for k = 1:Nperms
    comb2 = comb(randperm(length(comb)));
    xp = [xp, comb2(1:m)];
    yp = [yp, comb2(m+1:m+n)];
end



xfull=xstd*randn(m,Nperms);
yfull=randn(n,Nperms);



figure
subplot(121);hist(mean(xfull)-mean(yfull),Nperms/20)
subplot(122);hist(mean(xp)-mean(yp),Nperms/20)


%% Standard permutation test with bootstrapping up to even trial sizes
% Also produces incorrect distribution...

m=4;
n=100;
Nperms = 1500;
xstd=10;

x=xstd*randn(m,1);
y=randn(n,1);

%Bootstrap up x
x = x(floor(unifrnd(1,m+1-0.0000001,1,n)));

% Permutation test
comb = [x;y];
xp=[];
yp=[];
for k = 1:Nperms
    comb2 = comb(randperm(length(comb)));
    xp = [xp, comb2(1:n)];
    yp = [yp, comb2(n+1:n+n)];
end


xfull=xstd*randn(m,Nperms);
xfull2 = zeros(n,Nperms);
for k = 1:Nperms
    xfull2(:,k) = xfull(floor(unifrnd(1,m+1-0.0000001,1,n)),k);
end
yfull=randn(n,Nperms);


figure
subplot(121);hist(mean(xfull2)-mean(yfull),Nperms/20)
subplot(122);hist(mean(xp)-mean(yp),Nperms/20)


%% Pollard method - Compare histograms
m=4;
n=100;
Nperms = 150;
xstd=10;

x=xstd*randn(m,1);
y=randn(n,1);

% Bootstrap
xp = x(unidrnd(m,m*Nperms,1)); xp = reshape(xp,[m,Nperms]);
yp = y(unidrnd(n,n*Nperms,1)); yp = reshape(yp,[n,Nperms]);



zp = ( (mean(xp)-mean(yp)) - (mean(x) - mean(y))); 

figure; hist(zp,50); title(['Bootstrapped dist: std=' num2str(std(zp))]);

xfull=xstd*randn(m,Nperms);
yfull=randn(n,Nperms);

zptrue = mean(xfull)-mean(yfull);

figure; hist(zptrue,50);title(['Actual dist: std=' num2str(std(zptrue))]);



%% Pollard method - Test false positives

m=4;
n=20;
Nperms = 100;
xstd=1;
Ntests = 4000;

tests = false(1,Ntests);
% tic
for i = 1:Ntests

    x=xstd*randn(m,1);
    y=randn(n,1);

    % Bootstrap
    xp = x(unidrnd(m,m*Nperms,1)); xp = reshape(xp,[m,Nperms]);
    yp = y(unidrnd(n,n*Nperms,1)); yp = reshape(yp,[n,Nperms]);
  

    %zp = sqrt(n+m)*( (mean(xp)-mean(yp)) - (mean(x) - mean(y)));    % Seems too large
    zp = ( (mean(xp)-mean(yp)) - (mean(x) - mean(y))); % Instead try this...
    
    thresh = quantile(zp,0.9);
    tests(i) = (mean(x)-mean(y)) > thresh;
end
% toc

error_rate = sum(tests) / length(tests) *100
fprintf(['False positives under null distribution = ' num2str(error_rate) ' %. \n']);


