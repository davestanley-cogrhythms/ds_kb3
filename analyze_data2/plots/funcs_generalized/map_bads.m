

function [bads_out] = map_bads(mode1,mode2,mypairs1,mypairs2,bads)
    % sp is assumed to be in the form of mode1
    % we want to convert it to be of form of mode2

    data_type1 = mode2datatype(mode1);
    data_type2 = mode2datatype(mode2);
    
    if strcmp(data_type1,data_type2)                                % Same data types
        bads_out = bads;
    elseif strcmp(data_type1,'FFC') && strcmp(data_type2,'SFC')     % FFC to SFC
        good_pairs = mypairs1(~bads,:);
        good_unique_elects = unique(good_pairs(:));     % All electrodes considered by FFC pairs
        
        bads_out = true(1,length(mypairs2));
        for i = 1:length(mypairs2)
            bads_out(i) = sum( mypairs2(i,1) == good_unique_elects ) == 0;   % SFC element i is bad if it is not used in the list of good FFC electrodes
        end
        
    elseif strcmp(data_type1,'FFC') && strcmp(data_type2,'PSD')     % FFC to PSD
        good_pairs = mypairs1(~bads,:);
        good_unique_elects = unique(good_pairs(:));     % All electrodes considered by FFC pairs
        
        bads_out = true(1,length(mypairs2));
        for i = 1:length(mypairs2)
            bads_out(i) = sum( mypairs2(i,1) == good_unique_elects ) == 0;   % PSD element i is bad if it is not used in the list of good FFC electrodes
        end
        
    else
          error('Unknown combination of pair types');  
    end
    

end

