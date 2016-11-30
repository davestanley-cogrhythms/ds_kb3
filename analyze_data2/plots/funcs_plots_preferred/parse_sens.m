

function [sensp, bad_zeros] = parse_sens(parsestring,sens,sens_alt,spike_prefererred_curr,spike_preferred_alt,pls,cellnum)
    % s* - columns of sens
    % sa* - columns of sens_alt
    % p* - columns of pls
    
    silent_mode = 1;
    convert_spc_zeros_to_nans = 0;      % Setting this to one can introduce some errors in the parsing
                                        % Plus it is irrelevant when 0s are in the predictors, since setting data points to zero is the same as removing them
                                        %   anyways. That is changes in the fit coeficients (beta) won't affect the residuals (yi-yhati) associated
                                        %   with the zero terms, and therefore won't affect the MSE. Changes to the constant offset will affect the
                                        %   MSE, however.
                                        %
                                        % Note that this is not the case for the dependent variable (y), so should consider setting
                                        %   convert_spc_zeros_to_nans to 1 for y.
                                        % Overall, it's better to sue bad_zeros to remove unwanted datapoints
                                        
    mark_zeros_as_bad = 0;

    % Rename sens into shorthand variables ( sens(:,1) = s1; etc)
    for i = 1:size(sens,2)
        istr = num2str(i);
        evalstr = ['s' istr '= sens(:,' istr ');'];
        if ~silent_mode; fprintf ([evalstr '\n']); end
        eval(evalstr);
    end
    
    % Rename sens_alt into shorthand variables ( sens_alt(:,1) = sa1; etc)
    for i = 1:size(sens_alt,2)
        istr = num2str(i);
        evalstr = ['sa' istr '= sens_alt(:,' istr ');'];
        if ~silent_mode; fprintf ([evalstr '\n']); end
        eval(evalstr);
    end
    
    % Rename spike_prefererred_curr into shorthand variables ( spike_prefererred_curr(:,1) = spc1; etc)
%     if convert_spc_zeros_to_nans == 1
%         spike_prefererred_curr=double(spike_prefererred_curr);
%         spike_prefererred_curr(spike_prefererred_curr==0)=NaN;
%     end
    for i = 1:size(spike_prefererred_curr,2)
        istr = num2str(i);
        evalstr = ['spc' istr '= spike_prefererred_curr(:,' istr ');'];
        if ~silent_mode; fprintf ([evalstr '\n']); end
        eval(evalstr);
    end
    
    % Rename spike_preferred_alt into shorthand variables ( spike_prefererred_alt(:,1) = spa1; etc)
%     if convert_spc_zeros_to_nans == 1;
%         spike_preferred_alt=double(spike_preferred_alt);
%         spike_preferred_alt(spike_preferred_alt==0)=NaN;
%     end
    for i = 1:size(spike_preferred_alt,2)
        istr = num2str(i);
        evalstr = ['spa' istr '= spike_preferred_alt(:,' istr ');'];
        if ~silent_mode; fprintf ([evalstr '\n']); end
        eval(evalstr);
    end
    
    % Rename pls into shorthand variables ( pls(:,1) = p1; etc)
    for i = 1:size(pls,2)
        istr = num2str(i);
        evalstr = ['p' istr '= pls(:,' istr ');'];
        if ~silent_mode; fprintf ([evalstr '\n']); end
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
        if ~silent_mode; fprintf ([evalstr '\n']); end
        eval(evalstr);
    end

    sensp = zeros(size(sens,1),length(parsestring));
    for i = 1:length(parsestring)
        istr = num2str(i);
        evalstr = ['sensp(:,' istr ') = ' parsestring{i} ';'];
        if ~silent_mode; fprintf([evalstr '\n']);end
        eval(evalstr);
    end
    
    if mark_zeros_as_bad
        bad_zeros = any(sensp == 0,2);
    else
        bad_zeros = zeros(size(sensp,1),1);
    end
    
    bad_zeros = bad_zeros(:)';

end