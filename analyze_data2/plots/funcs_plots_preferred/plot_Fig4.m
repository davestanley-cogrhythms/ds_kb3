
function plot_Fig4(sr,sch)


    figure;
    set(gcf,'Position',[186 479 1292 477])
    sr_grouped = [sr(sch(:,1) == 1,[1:2]); sr(sch(:,2) == 1,[3:4])];
    sr_grouped = sort(sr_grouped,2);
    subplot(131); bar(mean(sr_grouped));
    hold on; errorbar(1:2, mean(sr_grouped),std(sr_grouped)/sqrt(size(sr_grouped,1)),'k.');
    title(['Preferred Category Scheme When Relevant']);
    ylabel('Firing Rate (normalized')
    [h p] = ttest(sr_grouped)

    sr_grouped = [sr(sch(:,1) == 1,[3:4]); sr(sch(:,2) == 1,[1:2])];        % Non preferred when relevant
    sr_grouped = sort(sr_grouped,2);
    subplot(132); bar(mean(sr_grouped));
    hold on; errorbar(1:2, mean(sr_grouped),std(sr_grouped)/sqrt(size(sr_grouped,1)),'k.');
    title(['Non-Preferred Category Scheme When Relevant']);
    xlabel('Morph Group')
    [h p] = ttest(sr_grouped)

    sr_grouped = [sr(sch(:,1) == 1,[5:6]); sr(sch(:,2) == 1,[7:8])];        % Non preferred when relevant
    sr_grouped = sort(sr_grouped,2);
    subplot(133); bar(mean(sr_grouped));
    hold on; errorbar(1:2, mean(sr_grouped),std(sr_grouped)/sqrt(size(sr_grouped,1)),'k.');
    title(['Preferred Category Scheme When Irrelevant']);
    [h p] = ttest(sr_grouped)


end