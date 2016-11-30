


function h = plot_scatterpls_and_groups(pls_stats,group,bad_any,indx,indy,group0)

%     ind1 = 1;
%     ind2 = 3;
    x = pls_stats(~bad_any,indx);
    y = pls_stats(~bad_any,indy);
    h = plott_fit(x,y,'k.');
    hold on; plot_scattergroups(x,y,[group],~bad_any);
    
    % Get legends
    
    if exist('group0','var')
        group_legend = Grp;
        group_legend(1).ctgs = indx;
        group_legend(2).ctgs = indy;
        group_legend = group_legend.query_legend(group0);
        xlabel(['FFC ' group_legend(1).legend]); ylabel(['FFC ' group_legend(2).legend]); % Might need fixing
    end
    
%     fprintf(['L1 Slope SchA Rel Vs Irrel= ' num2str(std(pls_stats(~bad_any,3))/std(pls_stats(~bad_any,1))) '\n']);
%     fprintf(['L1 Slope SchB Rel Vs Irrel= ' num2str(std(pls_stats(~bad_any,4))/std(pls_stats(~bad_any,2))) '\n']);
    
end