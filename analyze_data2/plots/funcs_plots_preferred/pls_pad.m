
function pls_new = pls_pad(pls,mode_subgroups,pad_mode)

    if nargin < 3
        pad_mode = 1;   % 1 - Pad infs after pls data, up to 11
                        % 2 - Pad infs between pls(:,:,1) and pls(:,:,2) up
                        %     until size 11. This maps 2-4 onto 9-11.
    end
    
    pad_value = 0;
    %pad_value = Inf;
    
    pls_new = pls;
    if size(pls,3) < 11 && mode_subgroups(2) == 0           % If pls is missing some ctgs (i.e. <11), fill it out with Infs
        switch pad_mode
            case 1

                    %Old, slow way
                    %temp = ones(size(pls(:,:,end))) * pad_value;
                    %pls_old = cat(3,pls(:,:,:),repmat(temp,[1,1,11-size(pls,3)]));
                    
                    %New faster way
                    sz = size(pls); sz(3) = 11;
                    pls_new = pad_value * ones(sz);
                    pls_new(:,:,1:size(pls,3)) = pls;
                    

            case 2
                    %Old, slow way
                    temp = ones(size(pls(:,:,end))) * pad_value;
                    pls_old = cat(3,pls(:,:,1),repmat(temp,[1,1,11-size(pls,3)]),pls(:,:,2:end));
                    

                    sz = size(pls); sz(3) = 11;
                    pls_new = pad_value * ones(sz);
                    pls_new(:,:,[1, end-size(pls,3)+2:end]) = pls;
        end
    end
    
end