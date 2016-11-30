
function [Tout Xout] = calc_ts_func_old(T,X,time_bin,fract_overlap,fract_maxgap,fhandle)
    
    plot_debug = 0;

    fract_shift = 1.0 - fract_overlap;
    shift = time_bin*fract_shift;
    
    t = T(:,1);
    bin_min = min(t):shift:(max(t)-time_bin);
    bin_max = bin_min + time_bin;

    dt = median(diff(t));
    
    Tout = zeros(length(bin_min),size(T,2));
%     Xout = zeros(length(bin_min),size(T,2));
    good_data = logical(zeros(length(bin_min),1));
    for i = 1:length(bin_min)
        index = (t >= bin_min(i) & (t < bin_max(i)));
        
        t_temp = t(index);
        fract_gap = 1.0 - length(t_temp)*dt ./ (time_bin-dt);
        
        if fract_gap <= fract_maxgap
            good_data(i) = 1;
            for j = 1:size(T,2)
                t2 = T(index,j);
                x2 = X(index,j);
                %Xout(i,j) = fhandle(x,y);
                Xout(i,j) = fhandle(t2,x2);
                Tout(i,j) = mean([bin_min(i) bin_max(i)]);
            end
        else
            good_data(i) = 0;
        end        
    end
    
    Xout = Xout(good_data,:);
    Tout = Tout(good_data,:);
    
    if plot_debug         
        figure;
        for i = 1:size(X,2)
            subplot(ceil((size(X,2)-1)/2),2,i);
            temp = max(abs([X(:,i)]));
            hold on; plot(T(:,i),X(:,i)/temp,'k.','MarkerSize',15);
            hold on; plot(Tout(:,i),Xout(:,i),'r.','LineWidth',2,'MarkerSize',25);
            title(['Freq band' num2str(i)]);
            add_stimseiz(ratN,'k','LineWidth',2);
        end
        
    end
end


