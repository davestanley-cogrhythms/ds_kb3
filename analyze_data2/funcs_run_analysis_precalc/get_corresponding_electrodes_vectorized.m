

function lfp2 = get_corresponding_electrodes_vectorized(units1,md,do_adjacent)
    Nunits = length(units1);
    lfp2 = zeros(1,Nunits);
    adj_missing = false(1,Nunits);
    if ~do_adjacent
        for i = 1:Nunits
            lfp2(i) = md.unit_LFP_pairs(2,units1(i));
        end

    else
        for i = 1:Nunits
            lfp2(i) = md.unit_LFP_pairs(2,units1(i));
            [lfp2(i), adj_missing(i)]= get_adjacent_electrode(md,units1(i));
        end


    end
    
end