

function plot_units_and_sfc(ns,sfc,filter,take_mean)
    if ~exist('take_mean','var'); take_mean = 1; end

    ns = ns(:,filter,:);
    sfc = sfc(:,filter,:);
    
    if ~take_mean    
        plott_ani_pairs(ns,@plot,sfc,@plot);
    else
        ns = squeeze(mean(ns,2));
        sfc = squeeze(mean(sfc,2));
        subplot(211); plot(ns);
        subplot(212); plot(sfc);
    end

end