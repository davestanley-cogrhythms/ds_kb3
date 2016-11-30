

% Plot reaction times comparison for congruent versus incongruent


    
%% Load md and conditions 

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

    

%% Plots for all files

clear rtc rtic
rt=[];
for i = 1:length(md)
    %% For loop
    statfun = @mean;
    
    mc = get_good_samples(match_congruent{i},md(i));
    mic = get_good_samples(match_incongruent{i},md(i));
    rts = md(i).match_rt;
    
    rt = [rt; statfun(rts(mc)) statfun(rts(mic))];
    rtc{i} = rts(mc);
    rtic{i} = rts(mic);
end

%% Plots
rtc=vertcat(rtc{:});
rtic=vertcat(rtic{:});


ind = 1:40; [p1,h] = signrank(rt(ind,1),rt(ind,2))
ind = 41:79; [p2,h] = signrank(rt(ind,1),rt(ind,2))

figure; plot(diff(rt,[],2) ./ min(rt,[],2)); xlabel('File Num'); ylabel('%change reaction time (msec)'); title(['monkeyL p = ' num2str(p1) ' monkey O p = ' num2str(p2)]);
figure; plot(rt); xlabel('File Num'); ylabel('ReactionTime (msec)'); title(['monkeyL p = ' num2str(p1) ' monkeyO p = ' num2str(p2)]);legend('Congruent','Incongruent')
figure; ind = 1:40; [n,b] = hist(rt(ind,1),10); bar(b,n); hold on; [n2,b2] = hist(rt(ind,2),10); bar(b2,n2,'r');  xlabel('ReactionTime (msec)'); ylabel('nTrials'); title(['p = ' num2str(p1)]); legend('Congruent','Incongruent')
figure; ind = 41:79; [n,b] = hist(rt(ind,1),10); bar(b,n); hold on; [n2,b2] = hist(rt(ind,2),10); bar(b2,n2,'r');  xlabel('ReactionTime (msec)'); ylabel('nTrials'); title(['p = ' num2str(p2)]); legend('Congruent','Incongruent')



figure; plot(rtc); hold on; plot(rtic,'r'); xlabel('Trial Num'); ylabel('ReactionTime (msec)'); title(['p = ' num2str(ranksum(rtc,rtic))]); legend('Congruent','Incongruent')
figure; [n,b] = hist(rtc,20); bar(b,n); hold on; [n2,b2] = hist(rtic,20); bar(b2,n2,'r');  xlabel('ReactionTime (msec)'); ylabel('nTrials'); title(['p = ' num2str(ranksum(rtc,rtic))]); legend('Congruent','Incongruent')



%% Formal plotting

clear group
group(1:4) = Grp;
i=0;
i=i+1;group(i).datastats = rt(1:40,1); group(i).legend = 'Congruent';
i=i+1;group(i).datastats = rt(1:40,2); group(i).legend = 'Incongruent';
i=i+1;group(i).datastats = rt(41:79,1); group(i).legend = 'Congruent';
i=i+1;group(i).datastats = rt(41:79,2); group(i).legend = 'Incongruent';

opts_PSC.hmask = blkdiag(true(2),true(2));
opts_PSC.paperfig_mode = 1;
fign; [h1, h, p, out.PSC] = plot_stats_custstruct(group,opts_PSC);


out.group = group;