function [bad_any,src,sra,pls,ns,ev,ssp] = remove_bads_data(bad_any,src,sra,pls,ns,ev,ssp)
    % Fill dataset values with nans for cells we know are bad!
    src(bad_any,:) = nan;
    sra(bad_any,:) = nan;
    pls(:,bad_any,:) = nan;
    ns(:,bad_any,:) = nan;
    ev(:,bad_any,:) = nan;
    ssp(:,bad_any,:,:) = nan;
    

end