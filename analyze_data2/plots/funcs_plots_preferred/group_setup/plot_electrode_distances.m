

function [hsp, group_with_distances] = plot_electrode_distances(group,md,mypairs,bad_any)
    % group_with_distances now contains electrode distances in datastats
    % field

    plot_raw_hists = 0;
    
    % Pull out electrode coordinates
    temp = [md.electrode_locations];
    AP = [temp.AP];
    ML = [temp.ML];
    
    % Also pull out file numbers and electrode numbers (used only for
    % debugging)
    elect_nums = [temp.electrode_num];  % Electrode numbers (not used; only for debugging)
    temp = 1:length(md);
    file_nums = arrayfunu(@(x,y) repmat(x,1,length(y.electrode_locations.electrode_num)),temp,md);
    fnm = horzcat(file_nums{:});        % File number index (not used; only for debugging)
    
    
    % Estimate distances for all good data
    ind1 = 1:length(mypairs);
    dists_all = sqrt( (AP(mypairs(ind1,1))-AP(mypairs(ind1,2))).^2 + (ML(mypairs(ind1,1))-ML(mypairs(ind1,2))).^2);
    dists_all = dists_all(~bad_any);
    
    
    
    % Calculate distances for each group
    for i = 1:length(group)
        % Generate array of electrode distances for each cell in group
        %mypairs_curr = group(i).metadata.mypairs;
        mypairs_curr = mypairs(any(group(i).cells,2),:);    % More traditional way of doing things
        dists_curr = sqrt( (AP(mypairs_curr(:,1))-AP(mypairs_curr(:,2))).^2 + (ML(mypairs_curr(:,1))-ML(mypairs_curr(:,2))).^2);
        group(i).datastats = dists_curr;
        
        % Estimate edges - we want number of edges based on dists_curr, but
        % edge range to span dists_all
        [N, edges] = histcounts(dists_curr);
        edges = linspace(min(dists_all),max(dists_all),10);
        
        % Calculate N for dists_curr
        [N, edges] = histcounts(dists_curr,edges);
        
        % Estimate N for all dists_all
        [Nall, edges] = histcounts(dists_all,edges);

        % Estimate edge centers
        len = length(edges)-1;
        temp = [1:len; 2:len+1]';
        bins = [edges(temp(:,1)); edges(temp(:,2))]';
        bins = mean(bins,2)';

        mypercent = N ./ Nall * 100;
        
        if ~plot_raw_hists
            figure; hsp = bar(bins,mypercent,'k');
        else
            
            figl; subplot(221); bar(bins,N); title('Dists curr');
            subplot(222); bar(bins,Nall); title('Dists all')
            subplot(223);
            hsp = bar(bins,mypercent,'k'); title('Percent electrodes')
            % Recalculate histogram for dists_all using auto-generated
            % edges
            subplot(224);
            [Nall, edges] = histcounts(dists_all);
            % Estimate edge centers
            len = length(edges)-1;
            temp = [1:len; 2:len+1]';
            bins = [edges(temp(:,1)); edges(temp(:,2))]';
            bins = mean(bins,2)';
            bar(bins,Nall);
            title('Dists all with edges auto-determined')
        end
        
    end
    
    group_with_distances = group;
end