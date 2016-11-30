

function spike_trace = get_spikes_n_range(units,range)
% Get spike bins in a range based (endings includsive)
% Allows for the possibility of multiple units occupying 
% a single index slot

% Inputs:
%     units - vector of unit indices
%     range - vector or matrix. If a vector it should be the form [A; B], where A and B are
%         the min and max ranges to search for spikes
%         if range is a matrix, then range is of the form [A(:)'; B(:)']
%         where A and B are vectors describing the max and min spikes to
%         pull out
% Outputs:
%     Produces a point process vector (intensity function) describing
%     the number of spikes at each index point

    plot_on = 0;
    N_spike_ranges = size(range,2);
    N_bins = range(2,1) - range(1,1) + 1;
    Nunits = length(units);
    
    if isvector(range)
        range=range(:);
    end

%     tic
%     bins = [min(min(range)):max(max(range))];
%     [spike_trace xout] = hist(units(:),bins(:));
%     
%     toc
    


%     % % Method 1 (for loop, slow)
%     tic
%     spike_trace = zeros((range(2,1) - range(1,1)+1),size(range,2));
%     for i = 1:size(range,2)
%         index = find((units >= range(1,i)) & (units <= range(2,i)));
%         units_in_range = units(index);
%         bins = range(1,i):range(2,i);
%         [spike_trace(:,i) xout]= hist(units_in_range,bins);
%         sp1=spike_trace;
%     end
%     toc
%     
%     
%     % % Method 2 - hist + histc combination
%     tic
%     range_matrix = 0:N_bins-1;
%     range_matrix = repmat(range_matrix(:),1,N_spike_ranges) + repmat(range(1,:),N_bins,1);
%     [spike_trace xout] = hist(units,range_matrix);
%     spike_trace = reshape(spike_trace,N_bins,N_spike_ranges);
%     
%     tops_edges = [range(1,:)-0.5; range(1,:)+0.5];
%     [tops] = histc(units,tops_edges(:));
%     tops = tops(1:2:end);
%     bottoms_edges = [range(2,:)-0.5; range(2,:)+0.5];
%     [bottoms] = histc(units,bottoms_edges(:));
%     bottoms = bottoms(1:2:end);
%     spike_trace(1,:) = tops(:)';
%     spike_trace(end,:) = bottoms(:)';
%     sp2 = spike_trace;
%     toc
%     
    % % Method 3 - histc combination - conceptually simpler and,
    % surprisingly, also faster, despite the fact that there are many more
    % bins! This seems to result from the fact that hist is slower than
    % histc.
%     tic
    
    range_matrix = bounds_to_indices(range);
    
    % Test
%     range_matrix2 = 0:N_bins-1;
%     range_matrix2 = repmat(range_matrix2(:),1,N_spike_ranges) + repmat(range(1,:),N_bins,1);
%     sum(sum(range_matrix ~= range_matrix2))
    
    
    
    myedges = [(range_matrix(:)-0.5)'; (range_matrix(:)+0.5)'];
    [spike_trace xout] = histc(units,myedges(:));
    spike_trace = spike_trace(1:2:end);     % Remove space inbetween our defined edges
    spike_trace = reshape(spike_trace,N_bins,N_spike_ranges);
    sp3 = spike_trace;
%     toc
    

%     
%     
%     test1 = [3 3 7 50 50 50 50 200];
% %     testrange = [2 3 4; 7 8 9; 200 201 202]';
%     testrange = [1.5 2.5 2.5 3.5 3.5 4.5];
%     [N xout] = histc(test1,testrange(:));


%    test1 = [51 51 52 53 56 58];
% %     testrange = [2 3 4; 7 8 9; 200 201 202]';
%     testrange = [50.5 51.5 55.5 56.5];
%     [N xout] = histc(test1,testrange(:));

    
%     
    if plot_on
        figure;
        n = 1:max(units);
        x = zeros(1,max(units));
        x(units) = 1;
        plot(n,x,'.');
        hold on;
        for i = 1:size(range,2)
            plot([range(1,i) range(2,i)],[0 0],'k','LineWidth',20);
        end
%         plot([range(1) range(1)], [0 1],'k','LineWidth',2);
%         plot([range(2) range(2)], [0 1],'k','LineWidth',2);
        
        plot(range_matrix,spike_trace,'ro');
    end
    

end





