function [sbar] = plots_regress(X0,names0,group,marker,opts)

    


    %% Setup dataset
%     close all
    plot_other_on = 0;
    plot_hist_responsevar = 0;      % Histogram of dependent variable (Y)
    plot_residuals_on = 0;
    plot_added = 0;
    plot_cook_pretrim = 0;
    plot_cook = 0;
    plot_bargraph = 1;
    plot_simulated_responses = 0;
    plot_simulated_residuals = 0;
    
    
    fitmode = 4;
        glmmode = {'gamma','reciprocal'};
%         glmmode = {'gamma','logit'};
%         glmmode = {'binomial','logit'};
        
        
    preproc_transform = 1;                  % 0-no transform; 1-zscore (DEFAULT); 2-log; 3-Drop half datapoints (bootstrap); 4-Add cell number predictor (1st Monkey only); 5-Add cell number predictor (both Monkeys)
    dependent_var_transform = 1;            % 0-use raw; 1-use preproc_transformed (usually z-scored; DEFAULT); 2- x100 for binomial glm; 3-fisher transform; 4-shifted fisher transform (seems to give superior Gaussian residuals)
    if fitmode == 5 || fitmode == 6
        if strcmp(glmmode{1},'gamma')       % Some overrides...
            dependent_var_transform = 0;    % Make use original ffc values, so dependent var is always +ve
        elseif strcmp(glmmode{1},'binomial')
            dependent_var_transform = 2;    % Scale FFC range up to 0-100, so that binomial distr doesn't default to bernoulli
        end
    end
    Nscale=100;
        
        
    remove_cook = 1;           % Removes high leverage points based on cook's distance (> 3x mean(cook))
    remove_leverage = 1;        % Removes high leverage points based on leverage (> 2*p/n)
    residtypeReg = 'Pearson';
    residtypeGLM = 'Deviance';
%          Options are:
%             GLM mdlglm.Residuals.Properties.VariableNames = {'Raw'    'LinearPredictor'    'Pearson'    'Anscombe'    'Deviance'}
%             Linear mdl.Residuals.Properties.VariableNames = {'Raw'    'Pearson'    'Studentized'    'Standardized''}


    % For robust fittting (Linear modeling) only.
    robust_ops.RobustWgtFun = 'bisquare'; % Default
    robust_ops.Tune = 4.685; % Defualt - approximately 95% as statistically efficient as the ordinary least-squares estimates; Decreasing the tuning constant increases the downweight assigned to large residuals;
    %robust_ops.Tune = 1.185;
    
    % Pull in options
    if nargin > 4
        % Opts is a structure for overriding default regression settings
        if ~isempty(opts)
            fprintf('Options struct specified - overriding default options \n');
            vars_pull(opts);
        end
    end
    
    switch residtypeGLM
        case 'Raw'; residtypeGLM2 = 'resid';
        case 'LinearPredictor'; residtypeGLM2 = 'nomatch';
        case 'Pearson'; residtypeGLM2 = 'residp';
        case 'Anscombe'; residtypeGLM2 = 'resida';
        case 'Deviance'; residtypeGLM2 = 'residd';
        otherwise
            residtypeGLM2 = 'nomatch';
    end
    
    switch fitmode
        case {1,2,3,4}
            residtype = residtypeReg;
        case 5
            residtype = residtypeGLM2;
        case 6
            residtype = residtypeGLM;
    end
    
    names=names0;
    X=X0;
    switch preproc_transform
        case 0              % Do nothing
        case 1
            X=zscore(X0); % Do zscore - this shouldn't matter
        case 2
            X=log(X0); % Log to downweight high amplitude data points (upweights low amplitude, though!)
    
        case 3
            N=size(X0,1);             % Keep half the datapoints at random
            keep=round(N/2);
            ind = 1:N;
            ind = ind(randperm(N));
            X=X0(ind(1:keep),:);

        case 4
            %Add variable for cell number for just the first monkey
            N=size(X0,1);
            cellnum1 = [1:182 zeros(1,N-182)]';   % Monkey 1
            %cellnum2 = [zeros(1,182) 1:N-182]';   % Monkey 2
            X = [X0 cellnum1]
            names = {names0{:} 'cellnum1'}
            
        case 5
            % Add single predictor for cell numbers of both monkeys
            N=size(X0,1);
            cellnum = [1:182 1:(N-182)]';
            X0 = [X0 cellnum]
            names = {names{:} 'cellnum'}

    end
    
    
    switch dependent_var_transform
        case 0
            % Use original
            X(:,end) = X0(:,end);
        case 1
            % Use preprocess transformed
        case 2
            % Setup data for binomial glm
            X(:,end) = Nscale*X0(:,end);
            %X(:,end+1) = 100;
        case 3
            % Use fisher transform
            X(:,end) = fisher_trans(X0(:,end));
        case 4
            % Use fisher transform
            X(:,end) = fisher_trans_shifted(X0(:,end));
    end
    
    if plot_hist_responsevar
        figure; Y=(X(:,end)); subplot(121); hist(Y,20); xlabel('Transformed coherence'); ylabel('N'); subplot(122); qqplot_with_stats(zscore(Y));
    end
    
    % Build dataset
    ds = dataset({X, names{:}});
    ds.Properties
    
    % Run the regression, at last
    to_exclude = [];
    [mdlc, resid] = dofit(ds,fitmode,glmmode,residtype,robust_ops,Nscale,to_exclude);
    
    % Remove outliers & redo fit
    if remove_cook; [larg] = id_larg_cook(mdlc); to_exclude = [to_exclude(:); larg];  end
    if remove_leverage; [larg] = id_larg_leverage(mdlc); to_exclude = [to_exclude(:); larg]; end
    to_exclude = [to_exclude(:); 33]; % Remove data point 33 because it's trouble
    [mdlc, resid] = dofit(ds,fitmode,glmmode,residtype,robust_ops,Nscale,to_exclude);
    
    % Old school code for calculating regression; built-in calculation of excludes
    %[mdlc2, resid] = dofit_old(ds,fitmode,glmmode,residtypeReg,residtypeGLM,residtypeGLM2,robust_ops,Nscale,remove_cook,remove_leverage);
    %diff([mdlc.Coefficients.pValue mdlc2.Coefficients.pValue],[],2)
    
    
    
    if plot_residuals_on
        plot_residuals(zscore(resid));
    end
    
    if plot_added
        %% Added variable plots
        cnames = mdlc.CoefficientNames;
        N=length(cnames);
        xp = (mdlc.Coefficients(:,4).pValue);
        %figure
        for i = 2:N
            if strcmp(cnames{i},'adj11') ; continue; end
            figure;
            plotAdded(mdlc,cnames{i});
            title(['N = ' num2str(size(ds,1)) ' p = ' num2str(xp(i))]);
        end
        %figure; plotAdded(mdl)
        %%
    end
    
    if plot_cook
        %% Leverage and Cook's Distance
        figure; plotDiagnostics(mdlc)
        figure; plotDiagnostics(mdlc,'cookd')
    end
    
    
    
    legend_arr=mdlc.CoefficientNames;
    x = double(mdlc.Coefficients(:,1).Estimate);
    xerr = double(mdlc.Coefficients(:,2).SE);
    xp = double(mdlc.Coefficients(:,4).pValue);
    h = zeros(size(xp));
    h(double(xp) < 0.05) = 1;
    h(double(xp) < 0.01) = 2;
    h(double(xp) < 0.001) = 3;

    ind = strcmp(legend_arr,'adj11') | strcmp(legend_arr,'(Intercept)');
    ind=~ind;
    
    if plot_bargraph
    %% Plot bargraph
        fign;
        [hbar herr] = barwitherr(xerr(ind),x(ind),'k');
        if any(cellfun(@length,legend_arr(ind)) > 3); xticklabel_rotate(1:length(legend_arr(ind)),90,legend_arr(ind), 'FontSize',get_fontlevel(1)*get_scaling_factor*0.5,'interpreter','none')  % If long legends, rotate
        else
            temp=legend_arr(ind);
            temp=strrep(temp,'_','\_');
            set(gca,'XTickLabel',temp);                                                     % Else plot normally
        end
    end
    % Pack into structure array, with each each array entry for a
    % different regressor
    sbar=struct('legend_arr',legend_arr,'x',num2cell(x'),'xerr',num2cell(xerr'),'h',num2cell(h'),'p',num2cell(xp'),'Formula',mdlc.Formula);
    
    
    if plot_simulated_responses
        %% Simulated responses
        Y = double(ds(:,end));
        if fitmode == 6 && strcmp(glmmode{1},'binomial')
            out = predict(mdlc,double(ds(:,1:end-1)),'BinomialSize',Nscale*ones(size(ds,1),1));
            out2 = random(mdlc,double(ds(:,1:end-1)),'BinomialSize',Nscale*ones(size(ds,1),1));
        else
            out = predict(mdlc,double(ds(:,1:end-1)));
            out2 = random(mdlc,double(ds(:,1:end-1)));
        end
        figure; plot([Y,out,out2]) 
        xlabel('data point'); ylabel('Reponse Value'); legend('Original','Predicted mean','Random simulation');
    end

    if plot_simulated_residuals
        %% Simulated residuals
        
        % Using the model we just fit, generate some new data
        out2 = [];
        in2 = [];
        Nrepeats=5;
        for i = 1:Nrepeats
            Xin = double(ds(:,1:end-1));
            if fitmode == 6 && strcmp(glmmode{1},'binomial')
                out2 = [out2; random(mdlc,Xin,'BinomialSize',Nscale*ones(size(ds,1),1))];
            else
                out2 = [out2; random(mdlc,Xin)];
            end
            in2 = [in2; Xin];
        end

        % Fit a 2nd model to this new data, and measure the residuals.
        % These are the expected residuals assuming the model is "true"
        ds2 = dataset({[in2, out2], names{:}});
        ds2.Properties
        to_exclude2 = [];
        [mdlsim, resid2] = dofit(ds2,fitmode,glmmode,residtype,robust_ops,Nscale,to_exclude2);
        
        
        % Plot residuals of simulated data = these should be close to gaussian
        plot_residuals(zscore_with_nans(mdlsim.Residuals.Raw));
        title('Residuals of simulated data');
        
        % Compare with raw data residuals
        figure;
        qqplot(zscore_with_nans(mdlc.Residuals.Raw),zscore_with_nans(mdlsim.Residuals.Raw));
        title('QQ plot residuals of real data vs simulated data');

    end
    
    %% Plots
    if plot_other_on
        
        
        

        
        %%
        figure; plotResiduals(mdlc,'histogram')         %DEfault
        %%
        figure; plotResiduals(mdlc,'probability')
        %%
        figure; plotResiduals(mdlc,'lagged')
        %%
        figure; 
        plotResiduals(mdlc,'fitted')
        %%
        figure; plotResiduals(mdlc,'caseorder') 
        %%
        figure; plotResiduals(mdlc,'symmetry') 
        %% Slice
        figure; 
        plotSlice(mdlc)
        %%
        figure; 
        plotEffects(mdlc)
        %%
        figure; 
        plotInteraction(mdlc,'s8','sa8','predictions')
        %%
        figure; 
        plotInteraction(mdlc,'SchARel','SchANR','predictions')
        
    end

end


function [larg] = id_larg_cook(model)

    
                 
    cookthresh = 3*mean(model.Diagnostics.CooksDistance);
    larg = find(model.Diagnostics.CooksDistance > cookthresh);

end



function [larg] = id_larg_leverage(model)

    leverage_coef = 2;      % Matlab default
    leverage_coef = 10;

    
    
    n = model.NumObservations;
    p = model.NumEstimatedCoefficients;

    
    leveragethresh = leverage_coef*p/n;
    larg = find(model.Diagnostics.Leverage > leveragethresh);
end



function [mdlc, resid] = dofit(ds,fitmode,glmmode,residtype,robust_ops,Nscale,to_exclude)
    
    switch fitmode
        case 1
            %% Fit model
            mdl = LinearModel.fit(ds);
            mdlc = mdl
            resid = mdlc.Residuals.(residtype);
        case 2
            %% Stepwise fit
%                     Removes non-essential predictors
            mdls = LinearModel.stepwise(ds);
            %mdls = LinearModel.stepwise(ds,'constant','upper','interactions');
            mdlc = mdls
            resid = mdlc.Residuals.(residtype);
        case 3
            %% Robust fitting
%                     Downweights outlier data points
%                     If tune is unspecified, the default value in the table is used.
%                     Default tuning constants give coefficient estimates that are
%                     approximately 95% as statistically efficient as the ordinary
%                     least-squares estimates, provided the response has a normal
%                     distribution with no outliers. Decreasing the tuning constant
%                     increases the downweight assigned to large residuals; increasing
%                     the tuning constant decreases the downweight assigned to large
%                     residuals.    
            mdlr = LinearModel.fit(ds,'RobustOpts',robust_ops);     % For methods, search Wiki for Iterative weighted least squares
            %mdlr = LinearModel.fit(ds,'plotmode ~ s1 + s2 + s1:s2','RobustOpts',robust_ops);
            
            mdlc = mdlr
            resid = mdlc.Residuals.(residtype);
            % resid = mdlc.Residuals.(residtype) .* mdlc.Robust.Weights;        % This doesn't work
            
        case 4
                        %% Robust fitting
%                     Downweights outlier data points
%                     If tune is unspecified, the default value in the table is used.
%                     Default tuning constants give coefficient estimates that are
%                     approximately 95% as statistically efficient as the ordinary
%                     least-squares estimates, provided the response has a normal
%                     distribution with no outliers. Decreasing the tuning constant
%                     increases the downweight assigned to large residuals; increasing
%                     the tuning constant decreases the downweight assigned to large
%                     residuals.    
            

            modelfun = @(ds,to_exclude) LinearModel.fit(ds,'linear','RobustOpts',robust_ops,'Exclude',to_exclude); % Default model function
            
            mdlrc = modelfun(ds,to_exclude);
            mdlc = mdlrc
            resid = mdlc.Residuals.(residtype);
          
        case 5
            %% GLM
            %mdlg = fitglm(ds,'ResponseVar','plotmode','interactions','distr','normal','link','identity')
            %[b,dev,stats] = glmfit(X(:,1:end-1),X(:,end),'normal','link','identity')
            
            [b,dev,stats] = glmfit(double(ds(:,1:end-1)),double(ds(:,end)),'gamma','Exclude',to_exclude); stats.p
            
            resid = stats.(residtype);
            
        case 6
            %% GLM using GLM class
%             mdlglm = fitglm(ds,'linear','distr','gamma')
%             mdlglm = fitglm(ds,'linear','distr','gamma','link','reciprocal')
%             mdlglm = fitglm(ds,'linear','distr','binomial','link','logit','BinomialSize',Nscale*ones(size(ds,1),1))
            if strcmp(glmmode{1},'binomial')
                mdlglm = fitglm(ds,'linear','distr',glmmode{1},'link',glmmode{2},'Exclude',to_exclude,'BinomialSize',Nscale*ones(size(ds,1),1));
            else
                mdlglm = fitglm(ds,'linear','distr',glmmode{1},'link',glmmode{2},'Exclude',to_exclude);
            end
            mdlc = mdlglm
            resid = mdlc.Residuals.(residtype);
            sum(isnan(resid)) 
            
            
            
    end
    
    if sum(isnan(resid)) > 0
        %warning('Nans in residuals....removing');
        %keyboard
        resid = resid(~isnan(resid));
        % Apparently fitglm doesn't insert nans into residual data, unlike fitlm.
        
    end
end


% % % % % % % % % % % % % % % % % % % % % % % OLD CODE, NOT USED % % % % % % % % % % % % % % % % % % % % % % %
% % function [mdlc, resid] = dofit_old(ds,fitmode,glmmode,residtype,residtypeGLM,residtypeGLM2,robust_ops,Nscale,remove_cook,remove_leverage)
% %     
% %     switch fitmode
% %         case 1
% %             %% Fit model
% %             mdl = LinearModel.fit(ds);
% %             mdlc = mdl
% %             resid = mdlc.Residuals.(residtype);
% %         case 2
% %             %% Stepwise fit
% % %                     Removes non-essential predictors
% %             mdls = LinearModel.stepwise(ds);
% %             %mdls = LinearModel.stepwise(ds,'constant','upper','interactions');
% %             mdlc = mdls
% %             resid = mdlc.Residuals.(residtype);
% %         case 3
% %             %% Robust fitting
% % %                     Downweights outlier data points
% % %                     If tune is unspecified, the default value in the table is used.
% % %                     Default tuning constants give coefficient estimates that are
% % %                     approximately 95% as statistically efficient as the ordinary
% % %                     least-squares estimates, provided the response has a normal
% % %                     distribution with no outliers. Decreasing the tuning constant
% % %                     increases the downweight assigned to large residuals; increasing
% % %                     the tuning constant decreases the downweight assigned to large
% % %                     residuals.    
% %             mdlr = LinearModel.fit(ds,'RobustOpts',robust_ops);     % For methods, search Wiki for Iterative weighted least squares
% %             %mdlr = LinearModel.fit(ds,'plotmode ~ s1 + s2 + s1:s2','RobustOpts',robust_ops);
% %             
% %             mdlc = mdlr
% %             resid = mdlc.Residuals.(residtype);
% %             % resid = mdlc.Residuals.(residtype) .* mdlc.Robust.Weights;        % This doesn't work
% %             
% %         case 4
% %                         %% Robust fitting
% % %                     Downweights outlier data points
% % %                     If tune is unspecified, the default value in the table is used.
% % %                     Default tuning constants give coefficient estimates that are
% % %                     approximately 95% as statistically efficient as the ordinary
% % %                     least-squares estimates, provided the response has a normal
% % %                     distribution with no outliers. Decreasing the tuning constant
% % %                     increases the downweight assigned to large residuals; increasing
% % %                     the tuning constant decreases the downweight assigned to large
% % %                     residuals.    
% %             
% %             plot_cook_pretrim=0;
% %             larg=[]; largc=[]; largl=[];
% %             modelfun = @(ds,larg) LinearModel.fit(ds,'linear','RobustOpts',robust_ops,'Exclude',larg); % Default model function
% %             
% %             if remove_cook; [largc, mdl_temp] = calc_larg_cook(modelfun,ds,[]); end
% %             if remove_leverage; [largl, mdl_temp] = calc_larg_leverage(modelfun,ds,[]); end
% %             
% %             larg = [larg; largc; largl; 33];      % Remove data point 33 because it's trouble
% %             
% %             if plot_cook_pretrim && (remove_leverage || remove_cook)
% %                 %% Plot Cook's pre removal
% %                 figure; plotDiagnostics(mdl_temp,'leverage')
% %                 figure; plotDiagnostics(mdl_temp,'cookd')
% %             end   
% %             
% %             
% %             larg(:)'
% %             mdlrc = modelfun(ds,larg);
% %             mdlc = mdlrc
% %             resid = mdlc.Residuals.(residtype);
% %           
% %         case 5
% %             %% GLM
% %             %mdlg = fitglm(ds,'ResponseVar','plotmode','interactions','distr','normal','link','identity')
% %             %[b,dev,stats] = glmfit(X(:,1:end-1),X(:,end),'normal','link','identity')
% %             
% %             [b,dev,stats] = glmfit(double(ds(:,1:end-1)),double(ds(:,end)),'gamma'); stats.p
% %             
% %             resid = stats.(residtypeGLM2);
% %             
% %         case 6
% %             %% GLM using GLM class
% % %             mdlglm = fitglm(ds,'linear','distr','gamma')
% % %             mdlglm = fitglm(ds,'linear','distr','gamma','link','reciprocal')
% % %             mdlglm = fitglm(ds,'linear','distr','binomial','link','logit','BinomialSize',Nscale*ones(size(ds,1),1))
% %             if strcmp(glmmode{1},'binomial')
% %                 mdlglm = fitglm(ds,'linear','distr',glmmode{1},'link',glmmode{2},'BinomialSize',Nscale*ones(size(ds,1),1));
% %             else
% %                 mdlglm = fitglm(ds,'linear','distr',glmmode{1},'link',glmmode{2});
% %             end
% %             mdlc = mdlglm
% %             resid = mdlc.Residuals.(residtypeGLM);
% %             
% %             
% %             
% %     end
% % end
% % 
% % 
% % function [larg, model] = calc_larg_cook(modelfun,ds,larg)
% % 
% %     model=modelfun(ds,larg);
% %                  
% %     cookthresh = 3*mean(model.Diagnostics.CooksDistance);
% %     larg = find(model.Diagnostics.CooksDistance > cookthresh);
% % 
% % end
% % 
% % 
% % 
% % function [larg, model] = calc_larg_leverage(modelfun,ds,larg)
% % 
% %     leverage_coef = 2;      % Matlab default
% %     leverage_coef = 10;
% % 
% %     model=modelfun(ds,larg);
% %     
% %     
% %     n = model.NumObservations;
% %     p = model.NumEstimatedCoefficients;
% % 
% %     
% %     leveragethresh = leverage_coef*p/n;
% %     larg = find(model.Diagnostics.Leverage > leveragethresh);
% % end
% % % % % % % % % % % % % % % % % % % % % % % END OLD CODE, NOT USED % % % % % % % % % % % % % % % % % % % % % % %


function plot_residuals(resid)
    %% Residuals Normal QQ Plot
    %[muhat,sighat]=normfit(resid)
    %PD = ProbDistUnivParam('normal', [muhat,sighat]);
    resid = resid(~isnan(resid));
    %resid = zscore(resid);
    PD = ProbDistUnivParam('normal', [0 1]); 
    figure; %qqplot(resid,PD);
    qqplot_with_stats(resid,PD);
    %figure; qqplot(resid);
    %figure; hist(resid, 50)
    %figure; qqplot(resid);
    [h, p, ksstat, cv] = kstest(resid);
%     warning('do AD test here too');
%     title(['kstest h=' num2str(h) ' p=' num2str(p) ]);

    if 1
        %% Compare CDF 
        figure;
%             subplot(121);
        [F,x] = ecdf(resid);
        plot(x,F);
        Ftheory = normcdf(x,0,1);
        hold on; plot(x,Ftheory,'r');
        dF = abs(Ftheory-F);
        [m, I] = max(dF);
        if Ftheory(I) > F(I)
            hold on; plot([x(I) x(I)],[Ftheory(I), Ftheory(I) - ksstat],'k','LineWidth',2);
        else
            hold on; plot([x(I) x(I)],[Ftheory(I), Ftheory(I) + ksstat],'k','LineWidth',2);
        end
        legend('Empirical CDF','Normal CDF','KS Stat');
        xlabel('Stat');ylabel('P');
        title(['KS Stat=' num2str(ksstat) ' p=' num2str(p)]);
%             subplot(122)
%             plot(x,dF);
%             hold on; plot([x(I) x(I)],[0 ksstat],'k','LineWidth',2)
    end
    %figure; plotResiduals(mdlc,'probability');
end

function out = zscore_with_nans(in)

    if ~isvector(in); error('Input must be 1D'); end
    
    out = zscore(in(~isnan(in)));

end