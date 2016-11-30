

function [varargout] = coherency_trials_downsample(varargin)
%     Usage: [varargout] = coherency_with_Cerr2_and_equal_trials(desired_Ntrials,baseline_subtract,fname,varargin)
%     Implements [varargout] = fname(varargin) while forcing the number of
%     trials to be desired_Ntrials. Coherency is biased based on the number
%     of trials present in the data. Hence, this command is used for
%     ensuring that the Ntrials bias is uniform across conditions.
% 
%     Like coherency_with_Cerr2, it allows you request return of  Cerr
%     regardless of whether you are actually running the Jackknife. Will
%     generate fake data otherwise (see coherency_with_Cerr2)
% 
%     Forms:
%         [varargout] = coherency_with_Cerr2_and_equal_trials(desired_Ntrials,baseline_subtract,fname,varargin)
%     Input:
%         fname - coherence function handle
%         desired_Ntrials - Number of trials to downsample to
%         baseline_subtract - To re-subtract baseline after downsampling trials
%         varargin - inputs to the coherence function (specified by fname)
%         varargout - outputs requested from fname. If fails, will contain a
%             final entry for Cerr
% 

    % Setup some parameters
    plot_debug = 0;
    downsample_mode = 2;    % Set downsample mode (1-3)
%         1 - Just drop random trials (wasteful)
%         2 & 3 - Generates multiple groupings of trials of size
%             desired_Ntrials. Then calculates the coherence of each
%             group and averages. Grouping methods described below:
%         2 - Randomly permutes trials and then draws groups of size 
%             desired_Ntrials until all trials used up. (No trials wasted
%             but maybe some randomness due to how trials are selected)
%         3 - bootstraps a large number of trial groupings of size
%             desired_Ntrials. There will be repition of trials.
    
    desired_Ntrials = varargin{1};
    baseline_subtract = varargin{2};
    
    % If desired_Ntrials is off
    if desired_Ntrials == 0
        [varargout{1:nargout}] = coherency_with_Cerr2(varargin{3:end});
    else    % Otherwise...
        
        % Setup data
        fname = varargin{3};
        data1 = varargin{4};
        data2 = varargin{5};
        if nargin >= 6; args = varargin(6:end);
        else args = {}; end
        
        % Make sure number of trials match and they are sufficient
        sz1 = size(data1,2);
        sz2 = size(data2,2);
        if sz1 ~= sz2; warning('Size mismatch between data1 and data2! Trouble is coming'); end
        if sz1 < desired_Ntrials;
            %warning (['Requested more trials than are actually in the data! Defaulting to coherency_with_Cerr2 (running with ' num2str(sz1)  ' trials.)']);
            [varargout{1:nargout}] = coherency_with_Cerr2(varargin{3:end});
%             for ii = 1:nargout
%                 varargout{ii} = varargout{ii}*Inf;  % Black ball conditions containing too few data points. Instead of 
%             end                                     % doing this though, can just detect in the plotting code and mark as bad
            return
        end
        
        % Generate collections of trials of size desired_Ntrials
        switch downsample_mode
            case 1      % - Just randomly pick some trials
                Niterations = 1;
                curr_trials = randi(sz1,desired_Ntrials,1);
            case 2      % - Randomly pick trials until all used up (most efficient option?)
                % Permutation mode 2 - get the number of iterations for 
                Niterations = floor(sz1 / desired_Ntrials);
                curr_trials = zeros(desired_Ntrials, Niterations);
                iter = randperm(sz1);
                for i = 1:Niterations
                    start = (1+(i-1)*desired_Ntrials);
                    stop = (i*desired_Ntrials);
                    curr_trials(:,i) = iter(start:stop);
                end
            case 3      % - Bootstrap many groupings of trials
                exclude_duplication = 1;
                Niterations = floor(sz1 / desired_Ntrials)*10;
                
                if exclude_duplication
                    curr_trials = zeros(desired_Ntrials, Niterations);
                    for j = 1:Niterations
                        temp = randperm(sz1);
                        curr_trials(:,j) = temp(1:desired_Ntrials);
                    end
                else % With duplication
                    curr_trials = randi(sz1,desired_Ntrials,Niterations);  
                end
        end
        %fprintf(['Ndesired= ' num2str(desired_Ntrials) '; Ntrials=' num2str(sz1) ' Niterations=' num2str(Niterations) '.\n']);
        
        % Loop through iterations, calculating coherences and storing in out_all
        clear out_all
        out_all = repmat({[]},1,nargout);
        for i = 1:Niterations
            
            data1_curr=data1(:,curr_trials(:,i),:);
            data2_curr=data2(:,curr_trials(:,i),:);
            
            data1_curr = do_baseline_subtract(data1_curr,baseline_subtract);
            data2_curr = do_baseline_subtract(data2_curr,baseline_subtract);
            
            [out_curr{1:nargout}] = coherency_with_Cerr2(fname,data1_curr,data2_curr,args{:});
            
            for j = 1:nargout
                out_all{j} = cat(4,out_all{j},out_curr{j});    % Stack the data along dimension 4.
            end
            
        end
        
        % Average outputs
        for j = 1:nargout
            if j~=2         % Not phi
                varargout{j} = mean(out_all{j},4);    % Average the data along dimension 4.
            elseif j == 2   % Phi - use circular statistics
                phi_all = out_all{j};
                sz=[size(phi_all,1),size(phi_all,2),size(phi_all,3),size(phi_all,4)];
                
                phi = zeros(sz(1),sz(2),sz(3));
                %tic
                for j2 = 1:sz(2)
                    for j3 = 1:sz(3)
                        phi(:,j2,j3) = circ_mean(squeeze(phi_all(:,j2,j3,:))');    % Average the data along dimension 4.
                    end
                end
                %toc
                
                varargout{j} = phi;
            end
            
            
        end
        
        if plot_debug
            Cave1 = varargout{1};
            Cave1_all = out_all{1};
            f1 = varargout{6};
            
            [varargout2{1:nargout}] = coherency_with_Cerr2(varargin{3:end});
            Cave2 = varargout2{1};
            f2 = varargout2{6};
            
            figure;
            hold on; plot(f2,Cave2,'k','LineWidth',2);
            hold on; plot(f1,Cave1,'b','LineWidth',2);
            hold on; plot(f1,squeeze(Cave1_all),'r:');
            hold on; plot(f2,Cave2,'k','LineWidth',2); %(Redraw)
            hold on; plot(f1,Cave1,'b','LineWidth',2); %(Redraw)
            legend(['Default, Ntrials=' num2str(sz1)],['Downsample & Average, Ntrials=' num2str(desired_Ntrials)],[num2str(Niterations) ' individual downsamples']);
            xlabel('f (Hz)'); ylabel('Coherence (0-1)');
            ylim([0 1]);
            
            
        end
        
        
    end
    
    


end


% 
% % Sample plotting code for comparing different downsampling modes
% 
% 
% 
% %%
% 
% cave_arr = [];
% 
% for i = 1:30
%     tic
%     [Cave,phi,S12ave,S1ave,S2ave,f,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencyc,currlfp1_ctg,currlfp2_ctg,params);
%     cave_arr = [cave_arr, Cave(:)];
%     toc
% end
% 
% %%
% cave_arr2 = [];
% 
% for i = 1:30
%     tic
%     [Cave,phi,S12ave,S1ave,S2ave,f,confC,phistd,Cerr]=coherency_trials_downsample(Ntrials_min,baseline_subtract,@coherencyc,currlfp1_ctg,currlfp2_ctg,params);
%     cave_arr2 = [cave_arr2, Cave(:)];
%     toc
% end
% 
% %% Plotting
% 
% figure; plot(cave_arr);
% figure; plot(cave_arr2);
% 
% 
% 
% figure; plott_matrix3D(cave_arr);
% figure; plott_matrix3D(cave_arr2);
