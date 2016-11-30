
function bootstrap_null_to_pvals(sfc_mode,curr_stage)


    %% Setup vars
    if ~exist('sfc_mode','var'); sfc_mode=2.201512; end
    if ~exist('curr_stage','var'); curr_stage=3; end

    %% First, get list of files to fix

    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix] = build_sfcmode(sfc_mode, mode_subgroups);
    outpath = fullfile(getpath('path_buffer_curr'),['test2mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(curr_stage)]);
    outpath2 = fullfile(outpath,'trimmed');
    mkdir_dave (outpath2);


    fn = dir(fullfile(outpath,'*.mat'));
    fn = {fn.name};


    %% Perform the operation
    parfor f = 1:length(fn)
        myconvert(outpath,fn,f,outpath2)
    end
    
    
    %% Rearrange things
    %keyboard
    system(['mkdir ' fullfile(outpath,'full')]);
    system(['mv ' fullfile(outpath,'*.mat') ' ' fullfile(outpath,'full')]);
    system(['mv ' fullfile(outpath2,'*.mat') ' ' fullfile(outpath)]);
    system(['rmdir ' fullfile(outpath2)]);
    

end


function [pvals, critval]= calcp(Cave1,Cave2,Cavep1,Cavep2, do_phi)

    dCave = Cave1 - Cave2;
    if do_phi; dCave = get_min_phi_diff(dCave); end
    dCave = abs(dCave);

    dCave_null = Cavep1 - Cavep2;
    if do_phi; dCave_null = get_min_phi_diff(dCave_null); end
        
    sz = size(dCave_null);
    %pvals = sum(dCave_null > abs(repmat(dCave,[1,1,sz(3)])),3) ./ (sz(3)*ones([sz(1:2)])); % One tailed
    pvals = sum( (dCave_null > abs(repmat(dCave,[1,1,sz(3)]))) | (dCave_null < -abs(repmat(dCave,[1,1,sz(3)]))) ,3) ./ (sz(3)*ones([sz(1:2)])); % Two tailed
    critval = prctile(dCave_null,95,3);
            
end
    

function z = get_min_phi_diff(z)
    
    ind=z>pi; z(ind)=z(ind)-2*pi; ind=z<-pi; z(ind)=z(ind)+2*pi;        % Bounds of maximum shift should be -pi to pi
    %z = abs(z);
    
end

function myconvert(outpath,fn,f,outpath2)
    clear sfc
    fprintf(['Loading file ' fullfile(outpath,fn{f}) '\n'])
    temp = load(fullfile(outpath,fn{f}));
    sfc = temp.sfc;

    for i = 1:length(sfc)

        s = sfc{i};
        
        % Calculate "difference" distribution, remove original dists to save space
        if isfield (s,'Cavep1')
            s.dCave_bs = s.Cavep1 - s.Cavep2;
            s = rmfield(s,{'Cavep1','Cavep2'});
        end
        if isfield(s,'phi1')
            s.dphi_bs = s.phip1 - s.phip2;
            s.dphi_bs = get_min_phi_diff(s.dphi_bs);
            s = rmfield(s,{'phip1','phip2'});
        end
        
        % Calculate statistics on Cave
        s.ks = calc_ks(s.Cave1,s.Cave2,s.dCave_bs,s.f,0);
        s.dCave_bs_sig = std(s.dCave_bs,[],3);
        s.zerotest = calc_zerotest(s.Cave1,s.Cave2,s.dCave_bs_sig,s.f);

        % Calculate statistics on phi
        if isfield(s,'phi1');
            s.phiks = calc_ks(s.phi1,s.phi2,s.dphi_bs,s.f,1);
            s.dphi_bs_sig = std(s.dphi_bs,[],3);
            s.phizerotest = calc_zerotest_phi(s.phi1,s.phi2,s.dphi_bs_sig,s.f);
            
        end
        
        if isfield (s,'dCave_bs')
            s = rmfield(s,{'dCave_bs'});
        end
        if isfield (s,'dphi_bs')
            s = rmfield(s,{'dphi_bs'});
        end
        
        

        sfc{i} = s;
    end
    save(fullfile(outpath2,fn{f}),'sfc');
end

function ks_struct = calc_ks(Cave1,Cave2,dCave_bs,f,isphi)
    %% calc_ks
    plot_on = 0;
    dCave = Cave1 - Cave2;
    if isphi; dCave = get_min_phi_diff(dCave); end % This variable dCave appears never to be used
    %dCave_bs = Cavep1 - Cavep2;
    
    sz = size(dCave_bs);
    h=false(sz(1:2));
    p=zeros(sz(1:2));
    ksstat=zeros(sz(1:2));
    for j = 1:sz(2)
        for i = 1:sz(1)
            freq = f(i);
            data = squeeze(dCave_bs(i,j,:));
            data = zscore(data);
            [h(i,j), p(i,j), ksstat(i,j), cv] = kstest(data);
        end
    end
    
    %ks_struct.h=h;
    ks_struct.p=p;
    ks_struct.ksstat=ksstat;
    
    if plot_on
        
        
        
        %% QQ and KS plots
        
        figure; imagesc(h); colorbar
        
        i=13;j=5;
        freq = f(i)
        data = squeeze(dCave_bs(i,j,:));
        data = zscore(data);
        PD = ProbDistUnivParam('normal', [0 1]);
        figure; qqplot(data,PD);
        %figure; qqplot(resid);
        %figure; hist(resid, 50)
        %figure; qqplot(resid);
        [h2, p2, ksstat2, cv] = kstest(data);
        title(['kstest h=' num2str(h2) ' p=' num2str(p2) ]);


        
        figure; %subplot(211);
        
        [F,x] = ecdf(data);
        plot(x,F);
        Ftheory = normcdf(x,0,1);
        hold on; plot(x,Ftheory,'r');
        dF = abs(Ftheory-F);
        [m, I] = max(dF);
        if Ftheory(I) > F(I)
            hold on; plot([x(I) x(I)],[Ftheory(I), Ftheory(I) - ksstat2],'k','LineWidth',2);
        else
            hold on; plot([x(I) x(I)],[Ftheory(I), Ftheory(I) + ksstat2],'k','LineWidth',2);
        end
        legend('Empirical CDF','Normal CDF','KS Stat');
        xlabel('Stat');ylabel('P');
        title(['KS Stat=' num2str(ksstat2) ' p=' num2str(p2)]);
        
%         subplot(212);
%         [n,xh] = hist(data,30);
%         bar(xh,n/sum(n));
%         hold on;
%         
    

    end
   

end


function zerotest = calc_zerotest(Cave1,Cave2,dCave_bs,f)
    % Percent deviation
    pd = 10;    % percent
    x = abs(Cave1 - Cave2);
    
    % Implement method to get p value
    sig = std(dCave_bs,[],3);
    mu = (Cave1 + Cave2) / 2;
    delta = mu * (pd/100);                % Reject null hypothesis for changes in Cave less than this value
    
    sz = size(x);
    p = 1 - normcdf(delta - x,zeros(sz),sig);
    zerotest.p = p;
    zerotest.pd=pd;
end


function zerotest = calc_zerotest_phi(Cave1,Cave2,dCave_bs,f)
    % Percent deviation
    angle_delta = 2*pi/10;    % allowed change in angle
    x = (Cave1 - Cave2);
    x = get_min_phi_diff(x);
    x = abs(x);

    
    % Implement method to get p value
    sig = std(dCave_bs,[],3);
    
    sz = size(x);
    p = 1 - normcdf(angle_delta - x,zeros(sz),sig);
    zerotest.p = p;
    zerotest.angle_delta=angle_delta;
end

