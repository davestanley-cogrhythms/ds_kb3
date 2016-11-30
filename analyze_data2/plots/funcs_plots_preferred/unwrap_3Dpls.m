

function pls_3D = unwrap_3Dpls(pls,f,f2)

    sz = size(pls);
    Nrows = length(f);
    Ncols = length(f2);
    
    pls_3D = reshape(pls,[Nrows,Ncols,sz(2:end)]);
    
    pls_3D = permute(pls_3D,[1,3,4,2]);
    
end