
function [wPEV, eta_sq] = extract_PEV(stage)
    % Calculates percent explained variation, unbiased estimator. See
    % supplement of Buschman, et al (2012) from Miller group.
    
    
    %% Calculate biased (eta_sq) and unbiased (omega squared AKA wPEV) percent explained variation
    plot_debug = 0;
    testing_mode = 0;

    path_prefunits = fullfile(getpath('path_buffer_mypreferred'));
    pref_prefix = 'preferred_';    
    prefname = [pref_prefix 'stage_' num2str(stage)];   

    pref_struct = load(fullfile(path_prefunits,prefname));
    
    % Spike rate, standard deviations
    sr = vertcat(pref_struct.spikerates{:});
    ss = vertcat(pref_struct.spikestd{:});
    Ntr = vertcat(pref_struct.sumctgsetli{:});
    Ncells = size(sr,1);
    N_ctgs = size(sr,2);
    N_pairings = floor((N_ctgs-1)/2); % Number of schemes
    
    % Testing mode
    if testing_mode
        Ncells = 5;
        sr=sr(1:Ncells,:);
        ss=ss(1:Ncells,:);
        Ntr=Ntr(1:Ncells,:);
        
    end
    
    % Sum squares
    SS = ss.^2 .* (Ntr - 1);
    

    
    % Hardcode some ctgs
    SchA_ctg = 9;
    SchB_ctg = 10;
    All_ctg = 11;
        
    % Calc sum of squares between groups
    SSbg = zeros(Ncells,N_pairings);
    i=0;
    i=i+1; gr1=1; gr2=2;  SSbg(:,i) = Ntr(:,gr1).*(sr(:,gr1)-sr(:,SchA_ctg)).^2 + Ntr(:,gr2).*(sr(:,gr2)-sr(:,SchA_ctg)).^2;      % Sch A Rel
    i=i+1; gr1=3; gr2=4;  SSbg(:,i) = Ntr(:,gr1).*(sr(:,gr1)-sr(:,SchB_ctg)).^2 + Ntr(:,gr2).*(sr(:,gr2)-sr(:,SchB_ctg)).^2;      % Sch B Rel
    i=i+1; gr1=5; gr2=6;  SSbg(:,i) = Ntr(:,gr1).*(sr(:,gr1)-sr(:,SchB_ctg)).^2 + Ntr(:,gr2).*(sr(:,gr2)-sr(:,SchB_ctg)).^2;      % Sch B NR
    i=i+1; gr1=7; gr2=8;  SSbg(:,i) = Ntr(:,gr1).*(sr(:,gr1)-sr(:,SchA_ctg)).^2 + Ntr(:,gr2).*(sr(:,gr2)-sr(:,SchA_ctg)).^2;      % Sch A NR
    i=i+1; gr1=9; gr2=10; SSbg(:,i) = Ntr(:,gr1).*(sr(:,gr1)-sr(:,All_ctg)).^2 + Ntr(:,gr2).*(sr(:,gr2)-sr(:,All_ctg)).^2;        % Variance in all trials explained by Scheme
    
    % Calc sum of squares total
    SStotal = zeros(Ncells,N_pairings);
    i=0;
    i=i+1; SStotal(:,i) = SS(:,SchA_ctg);  % All Sch A
    i=i+1; SStotal(:,i) = SS(:,SchB_ctg);  % All Sch B
    i=i+1; SStotal(:,i) = SS(:,SchB_ctg);  % All Sch B
    i=i+1; SStotal(:,i) = SS(:,SchA_ctg);  % All Sch A
    i=i+1; SStotal(:,i) = SS(:,All_ctg);   % All Ctg
    
    % Calculate MSE - their MSE is wrong - need to divide by error DOF. See http://davidmlane.com/hyperstat/B171693.html
    df_total = Ntr-1;
    df_groups = (2-1)*ones(size(Ntr));
    df_error = df_total - df_groups;
    
    
    MSE = zeros(Ncells,N_pairings);
    i=0;
    i=i+1; gr1=1; gr2=2;  MSE(:,i) = (SS(:,gr1)+SS(:,gr2)) ./ df_error(:,SchA_ctg);
    i=i+1; gr1=3; gr2=4;  MSE(:,i) = (SS(:,gr1)+SS(:,gr2)) ./ df_error(:,SchB_ctg);
    i=i+1; gr1=5; gr2=6;  MSE(:,i) = (SS(:,gr1)+SS(:,gr2)) ./ df_error(:,SchB_ctg);
    i=i+1; gr1=7; gr2=8;  MSE(:,i) = (SS(:,gr1)+SS(:,gr2)) ./ df_error(:,SchA_ctg);
    i=i+1; gr1=9; gr2=10; MSE(:,i) = (SS(:,gr1)+SS(:,gr2)) ./ df_error(:,All_ctg);
    
    eta_sq = SSbg ./ SStotal;                           % Biased PEV
    df = 2;
    wPEV = (SSbg - (df-1)*MSE) ./ (SStotal + MSE);      % Unbiased PEV
    
    
    %Try again for group 1
    top = Ntr(:,1).*(sr(:,1) - sr(:,9)).^2 + Ntr(:,2).*(sr(:,2) - sr(:,9)).^2;
    bottom = SS(:,9);
    MSE2 = SS(:,1) + SS(:,2);
    wPEV2 = (top - MSE2) ./ (bottom + MSE2);
    
    
    
    if plot_debug
        %% Debuggin plots
        figure;
        for i = 1:N_pairings
            clf; hist(wPEV(:,i),50)
            pause
        end
        
        
    end
    
end


function unused
%% Simple example code to implement example from davidmlane site
% Refs
% http://davidmlane.com/hyperstat/B86009.html
% http://davidmlane.com/hyperstat/B160638.html
% http://davidmlane.com/hyperstat/B171693.html


% Implement example
d = [3,5,3,5,2,2,4,4,2,1,3,2];
g1 = d(1:4);
g2 = d(5:8);
g3 = d(9:12);

df_total = length(d)-1;
df_groups = 3-1;
df_error = df_total-df_groups;

SST = sum((d - mean(d)).^2)
SSB = 4*(mean(g1) - mean(d)).^2 + 4*(mean(g2) - mean(d)).^2 + 4*(mean(g3) - mean(d)).^2
MSE = (sum((g1-mean(g1)).^2) + sum((g2-mean(g2)).^2) + sum((g3-mean(g3)).^2)) / df_error

eta_sq = SSB / SST
wPEV = (SSB - df_groups * MSE) ./ (SST + MSE)



end

