
% This function analyzes both sample and delay periods simultaneously! (as
% compared to run_preferred.m). Also saves my calculated preferred units.
function run_preferred3


    % % % Test Fries method on Kyle's simulated data. Then test Kyle's method
    % and compare
    % % % 
    format compact
    format long g

    run_setdefaultfig
    %if ~exist('stage','var'); stage = 3; end
    sfc_mode = 5;       % Mode for units - should always be 5
    recalc = 1; % Recalculating the total number of units. Probably don't need to do this.
    recalc2 = 1; % Recalculting the bulk of the code. Probably should do this unless just plotting.
    plot_on=1;
    [~, file_range] = estimate_filelist(5, 3);


    % Path to metadata
    metadata_path = fullfile(getpath('path_metadata'),'sansunits');
    file_list = get_filelist(metadata_path);
    file_list=file_list(file_range);
    Nfiles = length(file_list);

    % Detination for data output
    outpath = fullfile(getpath('path_buffer_mypreferred'));
    outname_prefix = 'preferred_';
    mkdir_dave(outpath);    


    % Load all metadata and also calculate the total number of units. This is
    % kinda slow, so inserted recalc clause.
    if recalc
        tot_Nunits=  0;
        for i = 1:Nfiles
            fprintf ('Counting units in file %d of %d, name %s \n',i,length(file_list),file_list{i});
            md{i} = load (fullfile(metadata_path,file_list{i}));     % Load metadata
            %No longer need to remove this since it's removed earlier in extract_metadata_merge.m
            %md{i} = rmfield(md{i},'unit_array'); % Remove unit fields which we don't need
            %md{i} = rmfield(md{i},'unit_ind');   % These take up a lot of space
            tot_Nunits = tot_Nunits + calc_Nunits(md{i});
        end
        save('temp.mat','md','tot_Nunits');
    else
        load ('temp.mat');
    end

    % Produce nicely formatted list of all units considered in this analysis.
    unit_list = md_to_unitlist(md,file_list);       % This is actually a list of file unit name pairs, but I am keeping it like this for legacy code


    % % Get preferred schemes for sample period
    sch_raw = [];
    sch_raw_irrel = [];
    sr = [];
    if ~get_iscromer
        stage = -2; save_pschemes (stage,file_list,tot_Nunits,outpath,outname_prefix);
        stage = -1; save_pschemes (stage,file_list,tot_Nunits,outpath,outname_prefix);
        stage = 0; save_pschemes (stage,file_list,tot_Nunits,outpath,outname_prefix);
    end
    stage = 2; save_pschemes (stage,file_list,tot_Nunits,outpath,outname_prefix);
    stage = 3; save_pschemes (stage,file_list,tot_Nunits,outpath,outname_prefix);
    stage = 6; save_pschemes (stage,file_list,tot_Nunits,outpath,outname_prefix);
    

end

function [schemes, spikerates, spikestd, funames, sum_ctgsetli, stagesir] = calc_unit_stats(filename,unit_path,tot_Nunits)
    % Schemes output is ctg1 ctg2 ctg1_nr ctg2_nr
    plot_on = 0;
    
    load(fullfile(unit_path,filename));
    Nunits = size(sfc{1}.spikerates,2);
    
    N_ctgs = length(sfc);
    N_pairings = floor((N_ctgs-1)/2); % Number of schemes
    
    schemes = zeros(Nunits,N_pairings);

    
    
    for j = 1:Nunits
        for i = 1:N_pairings
            alpha = 0.05 / tot_Nunits;
            alpha = 0.05; 
            %alpha = 0.05 / 4;           % Apparently Jefferson used 2 or 4 factor.
            
            %[schemes(j,i),p] = ttest2(sfc{i*2-1}.spikerates(:,j),sfc{i*2}.spikerates(:,j),'Alpha',alpha);
            [h,schemes(j,i)] = ttest2(sfc{i*2-1}.spikerates(:,j),sfc{i*2}.spikerates(:,j),'Alpha',alpha);
            if plot_on;
                group1 = i*2-1;
                group2 = i*2;
                [n1 x1] = hist( sfc{group1}.spikerates(:,j));
                [n2 x2] = hist( sfc{group2}.spikerates(:,j));
                figure;
                hold on; bar([x1; x2]',[n1; n2]');
                title([num2str(group1) ' vs ' num2str(group2) '. H=' num2str(schemes(j,i)) ' p=' num2str(p)]);
                legend('Group1','Group2');
                pause
                
            end
        end
    end
    
    spikerates = zeros(Nunits,N_ctgs);
    spikestd = zeros(Nunits,N_ctgs);
    for i = 1:N_ctgs
        spikerates(:,i) = mean(sfc{i}.spikerates);
        spikestd(:,i) = std(sfc{i}.spikerates);
    end
    
    sum_ctgsetli = repmat(sfc{1}.sum_ctgsetli,Nunits,1);
    stagesir = sfc{1}.stagesir;
    
    
    for i = 1:N_ctgs
        
    end
    
    
    for i = 1:N_ctgs
        Ntrials(:,i) = std(sfc{i}.spikerates);
    end
    
    if get_iscromer     % If we're working with Cromer data, add in empty columns for NR case ; these dummy variables will make subsequent analysis for roy and cromer easier
        dummy_sr = zeros(Nunits,4);
        dummy_sp = zeros(Nunits,2);
        spikerates = [spikerates(:,1:4) dummy_sr spikerates(:,5:end)];
        spikestd = [spikestd(:,1:4) dummy_sr spikestd(:,5:end)];
        schemes = [schemes(:,1:2) dummy_sp schemes(:,3:end)];
    end

    units_examined = sfc{1}.unit_names;
    units_examined=convert_unit_underscores(units_examined);
    funames = cellfunu(@(u) [filename ' ' u],units_examined);
end




function Nunits = calc_Nunits(md)
    
    Nunits = length(md.unit_names);

end

function plot_Fig4(sr,sch)


    figure;
    set(gcf,'Position',[186 479 1292 477])
    sr_grouped = [sr(sch(:,1) == 1,[1:2]); sr(sch(:,2) == 1,[3:4])];
    sr_grouped = sort(sr_grouped,2);
    subplot(131); bar(mean(sr_grouped));
    hold on; errorbar(1:2, mean(sr_grouped),std(sr_grouped)/sqrt(size(sr_grouped,1)),'k.');
    title(['Preferred Category Scheme When Relevant']);
    ylabel('Firing Rate (normalized')
    [h p] = ttest(sr_grouped)

    sr_grouped = [sr(sch(:,1) == 1,[3:4]); sr(sch(:,2) == 1,[1:2])];        % Non preferred when relevant
    sr_grouped = sort(sr_grouped,2);
    subplot(132); bar(mean(sr_grouped));
    hold on; errorbar(1:2, mean(sr_grouped),std(sr_grouped)/sqrt(size(sr_grouped,1)),'k.');
    title(['Non-Preferred Category Scheme When Relevant']);
    xlabel('Morph Group')
    [h p] = ttest(sr_grouped)

    sr_grouped = [sr(sch(:,1) == 1,[5:6]); sr(sch(:,2) == 1,[7:8])];        % Non preferred when relevant. This is correct (pair 1 with [5:6], etc)
    sr_grouped = sort(sr_grouped,2);
    subplot(133); bar(mean(sr_grouped));
    hold on; errorbar(1:2, mean(sr_grouped),std(sr_grouped)/sqrt(size(sr_grouped,1)),'k.');
    title(['Preferred Category Scheme When Irrelevant']);
    [h p] = ttest(sr_grouped)


end

function save_pschemes (stage,file_list,tot_Nunits,outpath,outname_prefix)
    % Get preferred schemes for a given stage
    % For output: Each cell is a file.
%     For pschemes, the 4 columns are:
%         SchemeA, SchemeB
%         
%     For spikerates, there are 9 columns. They are
%         8 columns for each of the 8 ctgs (ctgX and ctgX_nr)
%         1 column for total firing rate.
% 
%     For pschemes_irrel, the 4 columns are:
%         SchemeB irrelevant, SchemeA irrelevant
        
    
    sfc_mode = 5;
    Nfiles = length(file_list);

    unit_path = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode)],['stage_' num2str(stage)]);
    pschemes = {};
    spikerates = {};
    for i = 1:Nfiles 
        i
        fprintf ('Processing file %d of %d, name %s \n',i,length(file_list),file_list{i});
        [pschemes{i}, spikerates{i}, spikestd{i}, funames{i}, sumctgsetli{i}, stagesir{i}]=calc_unit_stats(file_list{i},unit_path,tot_Nunits);
    %     pschemes format is [a b c d], where:
    %         a=scheme A is significant;
    %         b=scheme B is significant
    %         c=scheme A differences are significant for when scheme B is present (irrelevant B)
    %         d=scheme B differences are significant for when scheme A is present (irrelevant A)
    
        pschemes{i}(isnan(pschemes{i})) = 0;    % If there are any nans in pschemes, it means the firing rate
                                                % for these cells is zero. In this case, we need to assume there
                                                % are no significant differences. (Will remove lowfiring cells later)
    end
    
    outname = [outname_prefix 'stage_' num2str(stage)];
    files_examined = file_list;
  
    save(fullfile(outpath,outname),'pschemes','spikerates','spikestd','files_examined','funames','sumctgsetli','stagesir');
    
end


function [sch_table labels] = calculate_cell_numbers(sch_raw)
    % Puts the raw preferred schemes in the format of Jefferson's table
    % Jefferson's labels:
%     A_sample_units_CRC
%     B_sample_units_CRC
%     A_delay_units_CRC
%     B_delay_units_CRC
% 
%     AB_sample_units_CRC
%     AB_delay_units_CRC
% 
%     A_sampledelay_units_CRC
%     B_sampledelay_units_CRC
%     AB_sampledelay_units_CRC
% 
%     AB_switch_dB_units_CRC: 0
%     AB_switch_sA_units_CRC: 0
%     AB_switch_sB_units_CRC: 0
%     AB_switch_dA_units_CRC: 0

    sch_table = zeros(size(sch_raw,1),13);
    
    % Returns TRUE if ALL of the columns are true and NO others
    i=0; clear labels 
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,1); labels{i} = 'A_sample';
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,2); labels{i} = 'B_sample';
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,3); labels{i} = 'A_delay';
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,4); labels{i} = 'B_delay';
    
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,[1 2]); labels{i} = 'AB_sample';
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,[3 4]); labels{i} = 'AB_delay';
    
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,[1 3]); labels{i} = 'A_sampledelay';
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,[2 4]); labels{i} = 'B_sampledelay';
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,[1 2 3 4]);  labels{i} = 'AB_sampledelay';
    
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,[1 2 3]); labels{i} = 'AB_sample_A_delay';
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,[1 2 4]); labels{i} = 'AB_sample_B_delay';
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,[1 3 4]); labels{i} = 'A_sample_AB_delay';
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,[2 3 4]); labels{i} = 'B_sample_AB_delay';
    
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,[1 4]); labels{i} = 'A_sample_B_delay';
    i=i+1; sch_table(:,i) = exclusive_true(sch_raw,[2 3]); labels{i} = 'B_sample_A_delay';

    
    % Should be a total of 15 rows in sch_table. This corresponds to 2^4-1 cases. The 16th case is nothing preferred
    % (the majority of units!)
    
    labels=labels(:);
end

function out = exclusive_true(in, columns)

    cols=zeros(1,4);
    cols(columns) = 1;
    columns_bar = find(~cols);
    
    out = and_across_rows(in(:,columns)) & ~or_across_rows(in(:,columns_bar));

end

function unit_list = md_to_unitlist(md,file_list)
    k = 0;
    for i = 1:length(md)
        for j = 1:length(md{i}.unit_names)
            k=k+1;
            unit_list{k} = [file_list{i} ' ' md{i}.unit_names{j}];
        end
    end
end


function royformat = generate_royformat(sch_table, labels, unit_list)

    royformat = [];
    Nfields = size(sch_table,2);
    for i = 1:Nfields
        index = find(sch_table(:,i));
        
        units_in_category = {};
        for j = 1:length(index)
           units_in_category{j} = unit_list{index(j)};
        end
        royformat = setfield(royformat,labels{i},units_in_category);
    end
end




