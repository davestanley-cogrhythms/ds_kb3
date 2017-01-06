

function pls = wrap_3Dpls(pls3D)
    % rearranges pls3D data so that frequency and timedata are both wrapped
    % along dimension 1. Originally, these data should be are stored along
    % dimensions 1 and 4, respectively. Undoes action of unwrap_3Dpls

    pls = permute(pls3D,[2,3,1,4]);
    pls = pls(:,:,:);
    pls = permute(pls,[3,1,2]);
    
end