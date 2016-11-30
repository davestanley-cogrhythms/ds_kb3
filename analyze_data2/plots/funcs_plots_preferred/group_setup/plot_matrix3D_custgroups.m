



function [h1, h2] = plott_matrix3D_custgroups(f,group1,group2,options,legends,fb,statsmode,stats_dozscore)

    mean3dmode=1;
    
        % Mode 1 = just mean across 3rd dim
        % Mode 2 = cat everything into second dim

    if size(group1,3) > 1
        if mean3dmode == 1
            group1=mean(group1,3);
            group2=mean(group2,3);
        else
            sz=size(group1); group1 = reshape(group1,[sz(1) sz(2)*sz(3)]); % (Can use reshape, it's probably better, works for sizes greater than 2!)
            sz=size(group2); group2 = reshape(group2,[sz(1) sz(2)*sz(3)]);
            %group1=cat(2,group1(:,:,1),group1(:,:,2));
            %group2=cat(2,group2(:,:,1),group2(:,:,2));
        end
    end

    if statsmode == 0
        h1 = plott_matrix3D(f,group1,options{:},'LineSpec',{'b-'});
        h2 = plott_matrix3D(f,group2,options{:},'LineSpec',{'g-'});
        legend([h1(1) h2(1)],[legends{1} ' N=' num2str(size(group1,2))],[legends{2} ' N=' num2str(size(group2,2))])
        xlabel('Freq (Hz)');
        ylabel('SFC - normalized');
    end
    
    index = f >= fb(1) & f < fb(2);
    X1=(group1);
    X2=(group2);
    if stats_dozscore
        X1=zscore(X1);
        X2=zscore(X2);
    end
    X1 = mean(X1(index,:,:),1);
    X1 = squeeze(X1); X1=X1(:);
    X2 = mean(X2(index,:,:),1);
    X2 = squeeze(X2); X2=X2(:);

    
    Xmu = [mean(X1) mean(X2)];
    Xste = [ste(X1,[],1) ste(X2,[],1)];
    [h p] = ttest2(X1,X2)
        
    if statsmode == 1


        [h1 h2 ] = bar_errsig(h,Xste,Xmu);
        ylabel(['SFC in Range ' num2str(fb(1)) '-' num2str(fb(2)) 'Hz range']);
    elseif statsmode == 2
        nbins = 20;
        subplot(121); hist(X1,nbins);
        subplot(122); hist(X2,nbins);
        h1=[];h2=[];
    end
    
    title(['p = ' num2str(p(1))])
    
    
end


function Y=ste(varargin)
    Y = std(varargin{:}) / sqrt(size(varargin{1},varargin{3}));
end
