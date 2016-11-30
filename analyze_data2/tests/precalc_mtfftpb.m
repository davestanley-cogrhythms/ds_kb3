

function precalc_mtfftpb(data1,data2,params)

[tapers,pad,Fs,fpass,err,trialave]=getparams(params);
N=check_consistency(data1,data2);
nfft=max(2^(nextpow2(N)+pad),N);
[f,findx]=getfgrid(Fs,nfft,fpass); 
tapers=dpsschk(tapers,N,Fs); % check tapers

%% Subtract in fourier domain

% Data 1
data=data1;
dbase = squeeze(mean(data,2));

J=mtfftc(data,tapers,nfft,Fs);
Jbase=mtfftc(dbase,tapers,nfft,Fs);

Jbase = repmat(Jbase,[1,1,size(J,3)]);
J1 = J-Jbase;

% Data 2
data=data2;
dbase = squeeze(mean(data,2));

figure; plott_matrix3D(data)
hold on; plot(dbase,'r');

J=mtfftpb(data,tapers,nfft);
Jbase=mtfftpb(dbase,tapers,nfft);
Jbase = repmat(Jbase,[1,1,size(J,3)]);
J2 = J-Jbase;

S12_fd=squeeze(mean(conj(J1).*J2,2));
S1_fd=squeeze(mean(conj(J1).*J1,2));
S2_fd=squeeze(mean(conj(J2).*J2,2));

% Take trial averages
S12_fd=squeeze(mean(S12_fd,2)); S1_fd=squeeze(mean(S1_fd,2)); S2_fd=squeeze(mean(S2_fd,2));



%% Subtract in time domain
% Data1
data = data1;
dbase = squeeze(mean(data,2));

dbase = repmat(dbase,[1,size(data,2)]);

data = data-dbase;
J1=mtfftc(data,tapers,nfft,Fs);

% Data2
data = data2;
dbase = squeeze(mean(data,2));

dbase = repmat(dbase,[1,size(data,2)]);

data = data-dbase;
J2=mtfftpb(data,tapers,nfft);


S12_td=squeeze(mean(conj(J1).*J2,2));
S1_td=squeeze(mean(conj(J1).*J1,2));
S2_td=squeeze(mean(conj(J2).*J2,2));


% Take trial averages
S12_td=squeeze(mean(S12_td,2)); S1_td=squeeze(mean(S1_td,2)); S2_td=squeeze(mean(S2_td,2));


%% Plot differences

figure;
plot(abs(S12_fd),'b')
hold on; plot(abs(S12_td),'ro')
legend('Freq domain calc','Time domain calc'); title('Abs(S12)');

figure;
plot(angle(S12_fd),'b')
hold on; plot(angle(S12_td),'ro')
legend('Freq domain calc','Time domain calc'); title('Phi(S12)');



figure;
plot(S1_fd,'b')
hold on; plot(S1_td,'ro')
legend('Freq domain calc','Time domain calc'); title('S1');


%
figure;
plot(S2_fd,'b')
hold on; plot(S2_td,'ro')
legend('Freq domain calc','Time domain calc'); title('S2');




