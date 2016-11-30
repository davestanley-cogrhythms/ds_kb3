

function group = get_group_overrides(fname_suffix, N_criteria, gr_submode)

    if ~exist('gr_submode','var'); gr_submode = []; end
    
    group = []; % Default output

    %N_criteria = 3;
    gr_all = Grp;
    gr_all.criteria = [2*ones(1,N_criteria)]; gr_all.criteria_alt = [2*ones(1,N_criteria)]; gr_all.ctgs = 1;


    % % Groups for reaction times; overrides others
    if ~isempty(strfind(fname_suffix,'RT')) && isempty(strfind(fname_suffix,'RT2'))
        clear group
        i=0;
        i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = 'Fast RT';
        i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = 'Slow RT';            

    end

    if ~isempty(strfind(fname_suffix,'RT2'))
        clear group
        i=0;
        if isempty(gr_submode); gr_submode = 1; end
        switch gr_submode
            case 2
                i=i+1; group(i) = gr_all; group(i).ctgs = [3]; group(i).legend = 'Incongruent Fast RT';
                i=i+1; group(i) = gr_all; group(i).ctgs = [4]; group(i).legend = 'Incongruent Slow RT';
            case 1
                i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = 'Congruent Fast RT';
                i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = 'Congruent Slow RT';
        end
    end

    if ~isempty(strfind(fname_suffix,'sw'))
        i=0; clear group
        i=i+1; group(i) = gr_all; group(i).ctgs = [5]; group(i).legend = 'Non-switch Trials';
        i=i+1; group(i) = gr_all; group(i).ctgs = [6]; group(i).legend = 'Switch Trials';
    end


    if ~isempty(strfind(fname_suffix,'bndry'))
        i=0; clear group
        if isempty(gr_submode); gr_submode = 1; end
        switch gr_submode
            case 3
                i=i+1; group(i) = gr_all; group(i).ctgs = [6]; group(i).legend = '50% morphs';
                i=i+1; group(i) = gr_all; group(i).ctgs = [5]; group(i).legend = '100% morphs';
            case 4
                i=i+1; group(i) = gr_all; group(i).ctgs = [8]; group(i).legend = '50% morphs irr';
                i=i+1; group(i) = gr_all; group(i).ctgs = [7]; group(i).legend = '100% morphs irr';
            case 1
                i=i+1; group(i) = gr_all; group(i).ctgs = [2]; group(i).legend = '60% morphs';
                i=i+1; group(i) = gr_all; group(i).ctgs = [1]; group(i).legend = '100% morphs';
            case 2
                i=i+1; group(i) = gr_all; group(i).ctgs = [4]; group(i).legend = '60% morphs irr';
                i=i+1; group(i) = gr_all; group(i).ctgs = [3]; group(i).legend = '100% morphs irr';
            case 5
                i=i+1; group(i) = gr_all; group(i).ctgs = [10]; group(i).legend = '50 or 60% morphs';
                i=i+1; group(i) = gr_all; group(i).ctgs = [9]; group(i).legend = '100% morphs';
            case 6
                i=i+1; group(i) = gr_all; group(i).ctgs = [12]; group(i).legend = '50 or 60% morphs irr';
                i=i+1; group(i) = gr_all; group(i).ctgs = [11]; group(i).legend = '100% morphs irr';
        end
    end

end