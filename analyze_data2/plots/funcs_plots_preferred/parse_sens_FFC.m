

function [sensp, bad_zeros] = parse_sens_FFC(parsestring,Lsens,Rsens,Lsens_alt,Rsens_alt,pls,cellnum,adj_stat)
    % sl* - columns of Lsens
    % sr* - columns of Rsens
    % sla* - columns of Lsens_alt
    % sra* - columns of Rsens_alt
    % p* - columns of pls
    %
    % cn1-2 - columns of cell numbers;
    % cn3-4 - monkey identifiers
    % cn5 - is percents switch
                                        
    mark_zeros_as_bad = 0;

    % Rename Lsens into shorthand variables ( Lsens(:,1) = sl1; etc)
    for i = 1:size(Lsens,2)
        istr = num2str(i);
        evalstr = ['sl' istr '= Lsens(:,' istr ');'];
        fprintf ([evalstr '\n']);
        eval(evalstr);
    end
        
    % Rename Rsens into shorthand variables ( Rsens(:,1) = sr1; etc)
    for i = 1:size(Rsens,2)
        istr = num2str(i);
        evalstr = ['sr' istr '= Rsens(:,' istr ');'];
        fprintf ([evalstr '\n']);
        eval(evalstr);
    end
    
    % Rename Lsens_alt into shorthand variables ( Lsens_alt(:,1) = sla1; etc)
    for i = 1:size(Lsens_alt,2)
        istr = num2str(i);
        evalstr = ['sla' istr '= Lsens_alt(:,' istr ');'];
        fprintf ([evalstr '\n']);
        eval(evalstr);
    end
    
    % Rename Rsens_alt into shorthand variables ( Rsens_alt(:,1) = sra1; etc)
    for i = 1:size(Rsens_alt,2)
        istr = num2str(i);
        evalstr = ['sra' istr '= Rsens_alt(:,' istr ');'];
        fprintf ([evalstr '\n']);
        eval(evalstr);
    end
    
    
    
    % Rename pls into shorthand variables ( pls(:,1) = p1; etc)
    for i = 1:size(pls,2)
        istr = num2str(i);
        evalstr = ['p' istr '= pls(:,' istr ');'];
        fprintf ([evalstr '\n']);
        eval(evalstr);
    end
    
    % Rename adj_stat into shorthand variables ( adj_stat(:,1) = adj1; etc)
    for i = 1:size(adj_stat,2)
        istr = num2str(i);
        evalstr = ['adj' istr '= adj_stat(:,' istr ');'];
        fprintf ([evalstr '\n']);
        eval(evalstr);
    end
    
    % Rename cellnum into shorthand variables ( cellnum(:,1) = cn1; etc)
    %  -- cn1 is the date of recording (days since beginning)
    %  -- cn2 is the date for Monkey L, with 0s for Monkey 0
    %  -- cn3 as cn2, with vice versa for Monkeys L/O
    %  -- cn4 is 1 for entries corresponding to Monkey L
    %  -- cn5 is percent_switch trials
    for i = 1:size(cellnum,2)
        istr = num2str(i);
        evalstr = ['cn' istr '= cellnum(:,' istr ');'];
        fprintf ([evalstr '\n']);
        eval(evalstr);
    end

    % Lastly, evaluate sensitivity matrix
    sensp = zeros(size(Lsens,1),length(parsestring));
    for i = 1:length(parsestring)
        istr = num2str(i);
        evalstr = ['sensp(:,' istr ') = ' parsestring{i} ';'];
        fprintf([evalstr '\n']);
        eval(evalstr);
    end
    
    if mark_zeros_as_bad
        bad_zeros = any(sensp == 0,2);
    else
        bad_zeros = zeros(size(sensp,1),1);
    end

    bad_zeros = bad_zeros(:)';
end