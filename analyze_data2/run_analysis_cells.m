

function run_analysis_cells(sfc_mode,curr_stage,pairsstart,pairsstop,i0,i0j0_split,do_precalc)

    % Convet i0 from string to number if necessary
    if isstr(i0)
        i0 = eval(i0);
    end
    
    % Decode sfc_mode
    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    sfc_mode_group = floor(sfc_mode); 
    mode_firstdecimal = mode_subgroups(1);
    mode_sixthdecimal = mode_subgroups(6);

    % Get full list of simulations
    load Nall.mat
    switch mode_firstdecimal
        case 0; Npairs_arr = Nunits;
        case 2; Npairs_arr = Nunits;
        case 3; Npairs_arr = Nunits_field;
        case 4; Npairs_arr = Nfield_field;
        case 5; Npairs_arr = Nunits_units;
        case 6; Npairs_arr = Nlfp;
        case 7; Npairs_arr = Nunits;
    end
    
    if i0j0_split == 0
        fi = pairsstart;        % Pairsstart actually codes for the file number when running in this mode.
        j0 = 1:Npairs_arr(fi);
        
        if do_precalc
            run_analysis_precalc(sfc_mode,curr_stage,fi,i0,j0,i0j0_split);
        else
            run_analysis(sfc_mode,curr_stage,fi,i0,j0,i0j0_split);
        end
        
%         if ~isGpuAvailable
%             run_analysis_precalc(sfc_mode,curr_stage,fi,i0,j0,i0j0_split);
%         else
%             run_analysis_precalc(sfc_mode,curr_stage,fi,i0,j0,i0j0_split);  % Use GPU code if available
%         end
        
    elseif i0j0_split == 1
        % Generate array of simulation indicies and corresponding file numbers
        % and pair numbers
        pairN = arrayfunu(@(x) 1:x, Npairs_arr);
        fileN = arrayfunu(@(x,y) x*ones(1,y), 1:length(Npairs_arr),Npairs_arr);
        pairN = horzcat(pairN{:});
        fileN = horzcat(fileN{:});
        
        % Run through desired set of pairs.
        for pairs_ind = pairsstart:pairsstop
            run_analysis(sfc_mode, curr_stage,fileN(pairs_ind),i0,pairN(pairs_ind),i0j0_split);
        end
    end

end


function OK = isGpuAvailable
    try
        d = gpuDevice;
        OK = d.SupportsDouble;
    catch
        OK = false;
    end
end
