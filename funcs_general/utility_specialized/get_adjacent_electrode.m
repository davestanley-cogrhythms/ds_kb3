
function [adj_electrode, adj_missing] = get_adjacent_electrode(md,curr_unit)

    use_backup = 1; % If the default adjacent electrode is missing, instaed use one of the pre-specified backup electrodes.

    curr_electrode = md.unit_LFP_pairs(2,curr_unit);
    
    lfp_names = md.lfp_names;
    N_electrodes = length(lfp_names);
    for i = 1:N_electrodes
        curr_name = lfp_names{i};
        elects(i) = str2num(strrep(curr_name,'PFC_Elect_',''));
    end
    
    curr_elec_num = elects(curr_electrode);
    adj_elec_num = swap_electrodes(curr_elec_num);
    
    adj_electrode = find(elects == adj_elec_num);
    
    % Check if the adjacent electrode was missing
    adj_missing = 0;
    if isempty (adj_electrode)
        adj_electrode = curr_electrode;
        adj_missing = 1;
        
        if use_backup
            adj_elec_num = get_backup_electrode(md, curr_unit, curr_elec_num);
            adj_electrode = find(elects == adj_elec_num);
            if isempty(adj_electrode);
                fprintf('Backup failed too!! Breaking. \n');
                pause
                return;
            end
        end
    end

end

function eout = swap_electrodes(e1)

    if mod(e1,2) == 0 e1 = e1 - 1;
    else e1 = e1 + 1;
    end
    
    eout = e1;

end


function adj_elec_num = get_backup_electrode(md, curr_unit, curr_elec_num);

    fname_md = strrep(md.recording_name,'-aligned.plx','');
    
    fname_md = [fname_md '.mat'];
    
    switch fname_md
        case 'L011107.mat'
            switch curr_unit
                case 5; adj_elec_num = 4;   % Currr elec PFC_Elect_7_unit_1
                case 6; adj_elec_num = 4;   % Currr elec PFC_Elect_7_unit_2
            end
        case 'L011707.mat'
            switch curr_unit
                case 4; adj_elec_num = 4;   % Currr elec PFC_Elect_7_unit_1
            end
        case 'L092906.mat'
            switch curr_unit
                case 7; adj_elec_num = 6;   % Currr elec PFC_Elect_7_unit_1
            end
        case 'O031805.mat'
            switch curr_unit
                case 4; adj_elec_num = 8;   % Currr elec PFC_Elect_6_unit_1
                case 7; adj_elec_num = 8;   % Currr elec PFC_Elect_16_unit_1
            end
            
        case 'O071105.mat'
            switch curr_unit
                case 4; adj_elec_num = 3;   % Currr elec PFC_Elect_5_unit_1
            end
    end
    
    if strcmp(fname_md(1:9),'ti0810161')
        switch curr_unit
            case 17; adj_elec_num = 14;   % Currr elec PFC_Elect_5_unit_1
            case 18; adj_elec_num = 14;   % Currr elec PFC_Elect_5_unit_1
        end
    end
    

end

