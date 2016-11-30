
function pairs = upper_triangle_pairs(Nelements)
        
        % Generate all possible pairs of triu matrix (off diagonal)
        [lfp1,lfp2] = meshgrid(1:Nelements,1:Nelements);
        lfp1 = tril(lfp1,-1);
        lfp2 = tril(lfp2,-1);
        pairs = [lfp1(:) lfp2(:)];
        pairs = pairs(lfp1 ~= 0,:);

end
