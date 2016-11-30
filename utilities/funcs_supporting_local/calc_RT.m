


function [tr_match,tr_nonmatch] = calc_RT(md)
    % Calculates the reaction times

    

end


% 
% function [sotr_match,sotr_nonmatch] = decode_match_nonmatch(md)
% 
%     % Set up codes - these are different for Cromer and Roy data!
% 
%     if ~get_iscromer
%         code_shift = 3;
%         code_sample_on = 27;
%         code_test = 29;
%         code_short = 4;
%         code_long = 30;
%     else
%         % Old method ... this one catches some bugs
%         %code_shift = 3;
%         %code_sample_on = 23;
%         %code_test = 25;
%         %code_short = 175;
%         %code_long = 26;
% 
%         code_shift = 4;
%         code_sample_on = 23;
%         code_test = 25;
%         code_short = 176;
%         code_long = 27;
%     end
%     
%     vars_pull(md)
%     ind27=find(all_encodes == code_sample_on);
% 
%     %Error trials
%     etrs = each_trial_response(sample_on_trials);
%     ind_err = etrs ~= 0;    % Trials for which the correct answer was not given (for whatever reason)
%                             % Not used!
% 
%     encodesp4 = all_encodes(ind27+code_shift);   % encode 3 steps down from sample on
%     sotr_match = encodesp4 == code_short;              % If this equals 4, it's a short (match) 
%     sotr_nonmatch = encodesp4 == code_long;             % If this equals 30, it's a long (non-match)
% 
% 
% end

