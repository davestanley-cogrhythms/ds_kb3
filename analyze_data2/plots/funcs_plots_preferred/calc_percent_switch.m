

%% Calc percent switch


function out = calc_percent_switch(md)

    clear Nsw Nct Ntrials* Ncells* percent_switch*
    for i = 1:79
       [tsw{i} tct{i}]= get_switch_trials(md(i));
       Nsw(i) = sum(tsw{i});
       Nct(i) = sum(tct{i});
       Ntrials_est(i) = Nsw(i)+Nct(i);
       Ntrials(i) = length(tsw{i});
       Ncells(i) = length(md(i).unit_names);
       Ncells2(i) = length(md(i).unit_LFP_pairs);
    end

    percent_switch = Nsw./(Ntrials) *100;


    %% Convert percent switch from one number per file to one number per cell
    percent_switch = Nsw./(Ntrials) *100;
    percent_switch=num2cell(percent_switch);

    for i = 1:length(Ncells)
        percent_switch{i} = repmat(percent_switch{i},1,Ncells(i));
    end
    percent_switch = horzcat(percent_switch{:});
    percent_switch = percent_switch(:);
    
    out.percent_switch=percent_switch;

end