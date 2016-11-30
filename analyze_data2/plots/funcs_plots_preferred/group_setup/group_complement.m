
function gout = group_complement(group,gr_all)

    %gr_all is only used to get the mapping from criteria to coords and vice versa
    allcoords = {gr_all(:).coords};
    allcoords = vertcat(allcoords{:});
    allcrit = {gr_all(:).criteria};
    allcrit = vertcat(allcrit{:});
    
    gout = group;
    for i = 1:length(group)
        gout(i).criteria(:,1:4) = gout(i).criteria(:,[2 1 4 3]);           % Swap criteria corresponding to A & B
        gout(i).coords = lookup_coords(gout(i).criteria,allcoords,allcrit);
        gout(i).ctgs = ctg_complement(gout(i).ctgs);
        
    end
    
    
    function ctgs_out = ctg_complement(ctgs)
        ctgs_out = ctgs;
        ind = ctgs == 1; ctgs_out(ind) = 3;
        ind = ctgs == 2; ctgs_out(ind) = 4;
        ind = ctgs == 3; ctgs_out(ind) = 1;
        ind = ctgs == 4; ctgs_out(ind) = 2;
        ind = ctgs == 5; ctgs_out(ind) = 7;
        ind = ctgs == 6; ctgs_out(ind) = 8;
        ind = ctgs == 7; ctgs_out(ind) = 5;
        ind = ctgs == 8; ctgs_out(ind) = 6;
    end
end



function coords = lookup_coords(crit,allcoords,allcrit)
    Ncrit = size(crit,1);
    coords = zeros(Ncrit, size(allcoords,2));
    for i = 1:Ncrit
        ind = find(all(repmat(crit(i,:),size(allcrit,1),1) == allcrit,2));  % Find index of criteria in allcrit matching crit(i,:)
        if isempty(ind); fprintf('Error: Complementary criteria not found. \n'); return; end
        coords(i,:) = allcoords(ind,:);
    end
end

