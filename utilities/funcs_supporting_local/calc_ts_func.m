
function [tout, Xout] = calc_ts_func(t,X,time_bin,fract_overlap,fract_maxgap,fhandle,cellularize_output)

%     [tout, Xout] = calc_ts_func(t,X,time_bin,fract_overlap,fract_maxgap,fhandle,cellularize_output)


    if ~exist('cellularize_output','var') cellularize_output=1;end
    
    plot_debug = 0;

    fract_shift = 1.0 - fract_overlap;
    shift = time_bin*fract_shift;
    
    if ~isvector(t); fprintf('t must be a vector'); tout=nan; Xout=nan; return; end     % Make into row vectors if necessary
    if isvector(X); X=X(:); end
    t=t(:);
    
    dt = mode(diff(t));
    
    bin_min = min(t):shift:(max(t)-time_bin);
    n_bin = round(time_bin/dt);
    bin_max = bin_min + time_bin;
    
    % Remove stupid rounding error which somehow creeps in
    bin_min = floor(bin_min*1e9)/1e9;
    bin_max = floor(bin_max*1e9)/1e9;

    % Test output of fhandle
    i=1;
    index = (t >= bin_min(i) & (t < bin_max(i)));
    [test_out]= apply_fhandle(t,X,fhandle,index);
    if ~cellularize_output    % If it is a scaler or a row vector, we can use a matrix. Else, use cell array
        outsize=size(test_out);
        outsize(1) = size(X,1);
        Xout = zeros(outsize);  % Preallocate matrix space to save memory access time.
    end
    
    tout = zeros(length(bin_min),1);
    good_data = false(length(bin_min),1);
    for i = 1:length(bin_min)
        %index = (t >= bin_min(i) & (t < bin_max(i)));
        index = find(t>=bin_min(i),1,'first');
        index = index:(index+n_bin-1);

        t_temp = t(index);
%         t_temp(1)
%         t_temp(end)
%         length(t_temp)
        fract_gap = 1.0 - length(t_temp)*dt ./ (time_bin-dt);
        
        if fract_gap <= fract_maxgap
            good_data(i) = 1;
            if ~cellularize_output
                Xout(i,:,:) = apply_fhandle(t,X,fhandle,index);
            else
                Xout{i,:,:} = apply_fhandle(t,X,fhandle,index);
            end
            tout(i) = mean([bin_min(i) bin_max(i)]);
        else
            good_data(i) = 0;
        end        
    end
    
    Xout = Xout(good_data,:,:);
    tout = tout(good_data);
    
    if plot_debug         
        figure;
        for i = 1:size(X,2)
            subplot(ceil((size(X,2)-1)/2),2,i);
            %temp = max(abs([X(:,i)]));
            hold on; plot(t,X(:,i),'k');
            hold on; plot(tout,Xout(:,i),'r');
            title(['Column' num2str(i)]);
        end
    end
end

function [Xout]= apply_fhandle(t,X,fhandle,index)
    Xout = fhandle(X(index,:,:));
    
%     if isnumeric(Xout) && isvector(Xout)
%         Xout = (Xout(:))';  % Make into a row so that it can be appended onto the output matrix appropriately.
%     end
end


