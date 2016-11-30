
function [src, ssc, pls, ns] = sortall_data(src,ssc,pls,ns,N_ctgs_base,N_ctgs_extra)
        % Making sort data function generalized...
    tic
    plot_debug = 0;

    N_ctgs = N_ctgs_base+N_ctgs_extra;

    N_pairs_base = round(N_ctgs_base/2);
    N_pairs_extra = round(N_ctgs_extra/2);
    N_pairs = N_pairs_base + N_pairs_extra;

    ctgs_pairs = 1:N_ctgs;
    ctgs_pairs = reshape(ctgs_pairs(:),2,round(N_ctgs/2));  % This matrix contains pairs of ctgs that will later be compared

    dom = {};
    for i = 1:N_pairs
        [~, dom{i}] = sort(src(:,[ctgs_pairs(:,i)]),2,'descend');
    end
    ctg1dom = dom{1};
    ctg3dom = dom{2};


    % % Test to make sure our sort_index_forloops function is working
    % X=src(:,[1:2]);
    % [Y inds] = sort(X,2);
    % Y2 = sort_index_forloops(X,2,inds);
    % Y3 = sort_index_forloops_old(X,2,inds);
    % sum(Y==Y2)
    % sum(Y==Y3)
    % length(Y)

    % % Sort based on dominant category

    % Plot ns groups post sorting
    if plot_debug
        %figure;plott_ani(ns(:,:,[1:2]))
        figure;
        subplotsq(3,1);plot(src(:,[1:2])); title('Post sorting - Firing rates')
        subplotsq(3,2);plot(squeeze(mean(pls(:,:,[1:2]),2))); title('SFC')
        subplotsq(3,3);plot(squeeze(mean(ns(:,:,[1:2]),2))); title('Traces')
    end


    % Units
    src(:,[1:2]) = sort_index_forloops(src(:,[1:2]),2,ctg1dom);     % Relevant Ctgs
    src(:,[3:4]) = sort_index_forloops(src(:,[3:4]),2,ctg3dom);
    if ~get_iscromer
        src(:,[5:6]) = sort_index_forloops(src(:,[5:6]),2,ctg1dom); % Non-relevant Ctgs (Roy only)
        src(:,[7:8]) = sort_index_forloops(src(:,[7:8]),2,ctg3dom);
    end
    for i = N_pairs_base+1:N_pairs
        curr_pairs = ctgs_pairs(:,i);                               % Extra Ctgs (schDom, Testing, match/nonmatch)
        src(:,curr_pairs) = sort_index_forloops(src(:,curr_pairs),2,dom{i});
    end
    
   % Units STD
    ssc(:,[1:2]) = sort_index_forloops(ssc(:,[1:2]),2,ctg1dom);     % Relevant Ctgs
    ssc(:,[3:4]) = sort_index_forloops(ssc(:,[3:4]),2,ctg3dom);
    if ~get_iscromer
        ssc(:,[5:6]) = sort_index_forloops(ssc(:,[5:6]),2,ctg1dom); % Non-relevant Ctgs (Roy only)
        ssc(:,[7:8]) = sort_index_forloops(ssc(:,[7:8]),2,ctg3dom);
    end
    for i = N_pairs_base+1:N_pairs
        curr_pairs = ctgs_pairs(:,i);                               % Extra Ctgs (schDom, Testing, match/nonmatch)
        ssc(:,curr_pairs) = sort_index_forloops(ssc(:,curr_pairs),2,dom{i});
    end

    % SFC traces
    pls = sortindex_all(pls,ctg1dom,ctg3dom,N_pairs_base,N_pairs,ctgs_pairs,dom);
    
    % Unit traces
    ns = sortindex_all(ns,ctg1dom,ctg3dom,N_pairs_base,N_pairs,ctgs_pairs,dom);



    % % Take frequency bands of SFC data for use in stats
    % index = f > freqband_stats(1) & f < freqband_stats(2);
    % if stats_dozscore; sfc_stats = zscore(sfc);
    % else sfc_stats = sfc;
    % end
    % sfc_stats=sfc_stats(index,:,:);
    % sfc_stats=mean(sfc_stats,1); sfc_stats = squeeze(sfc_stats);


    % Plot ns groups post sorting
    if plot_debug
        %figure;plott_ani(ns(:,:,[1:2]))
        figure;
        subplotsq(3,1);plot(src(:,[1:2])); title('Post sorting - Firing rates')
        subplotsq(3,2);plot(squeeze(mean(pls(:,:,[1:2]),2))); title('SFC')
        subplotsq(3,3);plot(squeeze(mean(ns(:,:,[1:2]),2))); title('Traces')
    end


    % Plot Miller Fig 4

    if plot_debug
        plot_Fig4(src,spc)
    end
    toc
end


function mat3D = sortindex_all(mat3D,ctg1dom,ctg3dom,N_pairs_base,N_pairs,ctgs_pairs,dom)

    mat3D(:,:,[1:2]) = sort_index_forloops(mat3D(:,:,[1:2]),3,permute(ctg1dom,[3 1 2]));
    mat3D(:,:,[3:4]) = sort_index_forloops(mat3D(:,:,[3:4]),3,permute(ctg3dom,[3 1 2]));
    if ~get_iscromer
        mat3D(:,:,[5:6]) = sort_index_forloops(mat3D(:,:,[5:6]),3,permute(ctg1dom,[3 1 2]));
        mat3D(:,:,[7:8]) = sort_index_forloops(mat3D(:,:,[7:8]),3,permute(ctg3dom,[3 1 2]));
    end
    for i = N_pairs_base+1:N_pairs
        curr_pairs = ctgs_pairs(:,i);                               % Extra Ctgs (schDom, Testing, match/nonmatch)
        mat3D(:,:,curr_pairs) = sort_index_forloops(mat3D(:,:,curr_pairs),3,permute(dom{i},[3 1 2]));
    end
end


