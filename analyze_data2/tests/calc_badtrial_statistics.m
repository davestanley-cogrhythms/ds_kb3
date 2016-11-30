

%% Calculate percent bads
bad_sum = cellfunu(@sum,bad_trials);
Ntrials = cellfunu(@(x) size(x,1),bad_trials);

percent_bad = cellfunu( @(x,y) x / y, bad_sum,Ntrials);

figure; plot(horzcat(percent_bad{:}));

%% Generate boxplot

pba = horzcat(percent_bad{:});  % Percent bads as single array

file_grouping = cellfunu( @(x,y) y*ones(1,length(x)),percent_bad,num2cell(1:length(percent_bad)));
file_grouping = horzcat(file_grouping{:});  % File numbers as single array

figure; boxplot(pba,file_grouping);





