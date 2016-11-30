

% As plots_reactiontimes, excepts compares RTs for boundary vs
% non-broundary trials

   

%% Load md and conditions 

run_setdefaultfig
addpath(genpath('./funcs_plots_preferred/'));


recalc_md=0;
reload_data = 0;

path_matfiles_store = getpath('path_matfiles_store');

% Create buffer variable
if ~exist('wrkspc_buffer','var'); wrkspc_buffer = struct; end

% Load Metadata
buff_fieldname = 'currmd';
buff_func = @() func_recalc_md(file_range);
wrkspc_buffer = buildstruct_buffered(wrkspc_buffer,path_matfiles_store, buff_fieldname, buff_func,recalc_md,reload_data);
md = wrkspc_buffer.(buff_fieldname).md;


condspath = getpath('path_conditions');
load(fullfile(condspath,'conditions_organized.mat'))



%% Import data if running as a function; define dependent parameters

running_as_script = ~is_calledby;
if ~running_as_script
    pop_workspace(false);
end

    

%% Calc for all files

clear rtb rtp
rt=[];
for i = 1:length(md)
    %% For loop
    statfun = @mean;
    
    [t60, t60_irr, t50, t50_irr, t80, t80_irr, t100, t100_irr] = get_boundary_trials(md(i));
%     tbndry = t60 | t50;
%     tbndry = t50;
    tbndry = t60;
    tpure = t100;
    rts = md(i).match_rt;
    
    rtb{i} = rts(tbndry);
    rtp{i} = rts(tpure);
    
    % Match trials only
    rtb{i} = rtb{i}(~isnan(rtb{i}));
    rtp{i} = rtp{i}(~isnan(rtp{i}));
    
    rt = [rt; statfun( rtb{i}) statfun( rtp{i})];
    
end

rtb=vertcat(rtb{:});
rtp=vertcat(rtp{:});

%% Plots
ind = 1:40; [p1,h] = signrank(rt(ind,1),rt(ind,2));
ind = 41:79; [p2,h] = signrank(rt(ind,1),rt(ind,2));

figure; plot(diff(rt,[],2) ./ min(rt,[],2)); xlabel('File Num'); ylabel('%change reaction time (msec)'); title(['monkeyL p = ' num2str(p1) ' monkey O p = ' num2str(p2)]);
figure; plot(rt); xlabel('File Num'); ylabel('ReactionTime (msec)'); title(['monkeyL p = ' num2str(p1) ' monkeyO p = ' num2str(p2)]);legend('Congruent','Incongruent')

warning('Ranksum is wrong here');
figure; ind = 1:40; [n,b] = hist(rt(ind,1),10); bar(b,n); hold on; [n2,b2] = hist(rt(ind,2),10); bar(b2,n2,'r');  xlabel('ReactionTime (msec)'); ylabel('nTrials'); title(['p = ' num2str(p1)]); legend('Boundary','Non-boundary')
figure; ind = 41:79; [n,b] = hist(rt(ind,1),10); bar(b,n); hold on; [n2,b2] = hist(rt(ind,2),10); bar(b2,n2,'r');  xlabel('ReactionTime (msec)'); ylabel('nTrials'); title(['p = ' num2str(p2)]); legend('Boundary','Non-boundary')


figure; plot(rtb); hold on; plot(rtp,'r'); xlabel('Trial Num'); ylabel('ReactionTime (msec)'); title(['p = ' num2str(ranksum(rtb,rtp))]); legend('Congruent','Incongruent')
figure; [n,b] = hist(rtb,20); bar(b,n); hold on; [n2,b2] = hist(rtp,20); bar(b2,n2,'r');  xlabel('ReactionTime (msec)'); ylabel('nTrials'); title(['p = ' num2str(ranksum(rtb,rtp))]); legend('Congruent','Incongruent')

%% Formal plotting

clear group
group(1:4) = Grp;
i=0;
i=i+1;group(i).datastats = rt(1:40,1); group(i).legend = '60% morphs';
i=i+1;group(i).datastats = rt(1:40,2); group(i).legend = '100% morphs';
i=i+1;group(i).datastats = rt(41:79,1); group(i).legend = '60% morphs';
i=i+1;group(i).datastats = rt(41:79,2); group(i).legend = '100% morphs';

opts_PSC.hmask = blkdiag(true(2),true(2));
opts_PSC.paperfig_mode = 1;
fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);


out.group = group;



