

function pls_3D = unwrap_3Dpls(pls,f,f2)
    % rearranges pls data so that frequency and timedata are stored along
    % dimensions 1 and 4, respectively. Originally, assumes pls has
    % frequency and time data wrapped along dimension 1.  Undoes action of wrap_3Dpls

    sz = size(pls);
    Nrows = length(f);
    Ncols = length(f2);
    
    pls_3D = reshape(pls,[Nrows,Ncols,sz(2:end)]);
    
    pls_3D = permute(pls_3D,[1,3,4,2]);
    
end