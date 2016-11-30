

function [sig_cells, pvalsf] = calc_significant_permutation_old_withplots(f,dCave,Cavep1,Cavep2,alpha,bad_any,freqval,do_bh_stepup,Cave)

    if ~exist('do_bh_stepup','var'); do_bh_stepup = 1; end      % Use BH setup up algorithm to enforce fixed false discovery rate, rather than just simply doing repeated t test
    plot_on = 0;
    plot_debug = 0;
    
    dCave = abs(dCave);
    dCave_null = Cavep1 - Cavep2;
        

    sz = size(dCave_null);
    %pvals = sum(dCave_null > abs(repmat(dCave,[1,1,sz(3)])),3) ./ (sz(3)*ones([sz(1:2)])); % One tailed
    pvals = sum( (dCave_null > abs(repmat(dCave,[1,1,sz(3)]))) | (dCave_null < -abs(repmat(dCave,[1,1,sz(3)]))) ,3) ./ (sz(3)*ones([sz(1:2)])); % Two tailed
    critval = prctile(dCave_null,95,3);
    
    
    
    % Look only at a specific frequency
    index = find(f >= freqval, 1, 'first');
    dCf = squeeze(dCave(index,:));
    dCf_null = squeeze(dCave_null(index,:,:));
    pvalsf = squeeze(pvals(index,:));
    critvalf = squeeze(critval(index,:));

    if do_bh_stepup
        
        [sig_cells, p0] = bh_step1D(pvalsf,alpha,bad_any);    % My code
        %sig_cells2 = fdr_bh(pvals,alpha);  % Mathworks code - they seem to give the same result
        
        
    else
        sig_cells = pvalsf < alpha;
        p0 = alpha;
    end
    
    
    
    %% Pull out cells for plotting
    plotCells = find(sig_cells);
    plotCells = plotCells(1:min(10,end));
    
    if plot_on
        %% Plot distributions of permuted spectra and null distribution
        figure; h1=plott_matrix3D(f,Cavep1(:,plotCells,:),'active_dim',3,'do_shift',0.01);
        hold on; h2=plott_matrix3D(f,Cavep2(:,plotCells,:),'active_dim',3,'do_shift',0.01);
        legend([h1(1),h2(1)],'Perm Ctg1','Perm Ctg2')
        title('Perm Ctg1 vs Perm Ctg2')

        
        figure; h1 = plott_matrix3D(f,dCave_null(:,plotCells,1:200),'active_dim',2,'do_shift',0.1,'do_mean',0,'showErrorbars',0)
        hold on; h2 = plott_matrix3D(f,dCave(:,plotCells),'active_dim',3,'do_shift',0.1,'LineSpec',{'LineWidth',2,'color','k'})
        legend([h1{1}(1),h2(1)],'Null','Stat')
        title('Null vs Stat');
    end    
    
    
    if plot_on
        %% Plot regions of spectra exceeding significance
        
        dCave_pass = dCave; dCave_pass(pvals > 0.05) = NaN;
        
        figure; h1=plott_matrix3D(f,dCave(:,plotCells),'active_dim',2,'do_shift',0.1,'do_mean',0,'showErrorbars',0)
        hold on;h2=plott_matrix3D(f,critval(:,plotCells),'active_dim',2,'do_shift',0.1,'do_mean',0,'showErrorbars',0,'LineSpec',{'LineWidth',2,'color','k'})
        hold on; h3=plott_matrix3D(f,dCave_pass(:,plotCells),'active_dim',2,'do_shift',0.1,'do_mean',0,'showErrorbars',0,'LineSpec',{'LineWidth',2,'color','r'})
        title('Stat vs crit value')
        legend([h1{1},h2{1},h3{1}],'Stat','Stat Crit value','Stat passes test')
        
        
    end
    
    
    if plot_on
        %% Plot histograms of null distribution versus data
        figure
        hsp=subplot_grid(length(plotCells))
        for i = 1:length(plotCells)
            hsp.set_gca(i)
            hist(dCf_null(plotCells(i),:),50);
            hold on;
            plot([dCf(plotCells(i)),dCf(plotCells(i))],[0 10],'k','LineWidth',5)
        end
    end
    
    
    if plot_on
        %% Plot significantly different spectra
        figure;
        hsp = subplot_grid(length(plotCells));
        Nnull = size(dCave_null,3);
        for i = 1:length(plotCells)
            hsp.set_gca(i)
            plot(f,squeeze(dCave_null(:,plotCells(i),1:Nnull)),'k:');
            hold on; plot(f,squeeze(Cave(:,plotCells(i),[9 10])),'LineWidth',2);
        end
        
        %figure; h1 = plott_matrix3D(f,dCave_null(:,plotCells,1:200),'active_dim',2,'do_shift',0.1,'do_mean',0,'showErrorbars',0)
        
        
    end
    
end


%%