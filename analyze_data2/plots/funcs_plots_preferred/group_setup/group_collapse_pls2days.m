

function group = group_collapse_pls2days(group)

    for i = 1:length(group)
        pls = group(i).data;
        pls_stats = group(i).datastats;
        funames1D = group(i).funames;
        
        sz=size(pls);
        ids = cellfun(@fname_to_filedate,funames1D);    % ID's
        idsu = unique(ids);                             % ID's unique
        Ndays=length(idsu);
        pls_new    =    zeros([size(pls,1),Ndays,size(pls,3)]);
        pls_stats_new = zeros([1,Ndays]);
        
        for j = 1:Ndays
            ind_today=ids==idsu(j);
            pls_new(:,j,:) = mean(pls(:,ind_today ,:),2);
            pls_stats_new(j) = mean(pls_stats(1,ind_today,:,:),2);
            funames_new{j} = funames1D(ind_today(1));
        end
        
        group(i).data = pls_new;
        group(i).datastats = pls_stats_new;
        group(i).numcells = Ndays;
        group(i).funames = funames_new;
        
    end


end
