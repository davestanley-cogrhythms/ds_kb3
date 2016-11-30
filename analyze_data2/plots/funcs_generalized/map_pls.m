

function [pls_mapped, bad_mapped] = map_pls(mode_src,mode_dest,pls,bad,mypairs_src,mypairs_dest)
    % sp is assumed to be in the form of mode1
    % we want to convert it to be of form of mode2
    
    plot_debug = 0;

    %data_type1 = mode2datatype(mode1);
    %data_type2 = mode2datatype(mode2);
    
    [pairs_type1, pairs_mode1 ]= mode2pairstype(mode_src);
    [pairs_type2, pairs_mode2 ]= mode2pairstype(mode_dest);
    
    fprintf(['Converting ' pairs_type1 ' to ' pairs_type2 '\n']);
    
    if pairs_mode1 == pairs_mode2
        %% Same data types
        pls_mapped = pls;
        bad_mapped = bad;
        
    elseif ( pairs_mode1 == 6 ) && pairs_mode2 == 4     % PSD to FFC
        pls_mapped = sqrt(pls(:,mypairs_dest(:,1),:) .* pls(:,mypairs_dest(:,2),:));    % Sqrt product of the two.
        bad_mapped = bad(mypairs_dest(:,1)) | bad(mypairs_dest(:,2));
    else
        %%
        error('Unknown combination of pair types');  
    end
end

