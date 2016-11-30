
function [sens, sensnames, ms_funames] = calc_sens(src,ssc,curr_stage,use_morph_sens,sens_std_normalize)

if use_morph_sens
    [sens,ms_funames] = loadall_morph_sens(curr_stage);
else
    sens = [src(:,1)-src(:,2), src(:,3)-src(:,4), src(:,5)-src(:,6), src(:,7)-src(:,8)]; % Estimate sensitivities by differences in firing rates
    %sens = [src(:,1)./src(:,2), src(:,3)./src(:,4), src(:,5)./src(:,6), src(:,7)./src(:,8)]; % Estimate sensitivities by differences in firing rates
end

% Add Test conditions
%sens = [sens, src(:,9)-src(:,10), src(:,11)-src(:,12)];

% Add Scheme sensitivity
sens = [sens, src(:,9)-src(:,10)];

% Take abs - might want to move this up!
sens = abs(sens);

if sens_std_normalize
    sens_stds = [ssc(:,1).^2+ssc(:,2).^2, ssc(:,3).^2+ssc(:,4).^2, ssc(:,5).^2+ssc(:,6).^2, ssc(:,7).^2+ssc(:,8).^2];
    %sens_stds = [sens_stds, ssc(:,9).^2+ssc(:,10).^2, ssc(:,11).^2+ssc(:,12).^2];
    sens_stds = [sens_stds, ssc(:,9).^2+ssc(:,10).^2];
    sens_stds = sqrt(sens_stds);
    
    sens = sens ./ sens_stds;
end

sens_main = sens(:,1:4);      % Sample image Sch A B Rel / NRel
sens_other = sens(:,5:end);   % Test image Sch A B Rel


senscirc = [sqrt((sens(:,1)+sens(:,2)).^2) sqrt((sens(:,3)+sens(:,4)).^2)];     % Lets each sensitivity measure be the radius of a circle passing through sens(:,1) and sens(:,2)
senscircdiff = senscirc(:,1) - senscirc(:,2);
sensprod = sqrt(sens(:,1).*sens(:,2));
sens = [sens_main, senscirc, senscircdiff sensprod src(:,end), sens_other];
sensnames = {'SchAR' 'SchBR' 'SchBnR' 'SchAnR' 'Circ_AB_R' 'Circ_AB_nR' 'Circ_R-nR', 'SensProd','FR','SchAB'};

end

