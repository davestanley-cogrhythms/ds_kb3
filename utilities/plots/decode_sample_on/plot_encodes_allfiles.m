

% This code is used to verify that the convert_to_seconds part of
% extract_metadata_merge is working correctly.


clear
close all

addpath(genpath('../../funcs_supporting_local/'));

path_md = fullfile(getpath('path_metadata'),'sansunits');
%path_md = fullfile(getpath('path_metadata'));
%path_md = fullfile(getpath('path_roydata_orig'),'metadata');
%path_md = getpath('path_roy_newctgX');

file_list = get_filelist(path_md);

md = cellfun(@load,fullfile(path_md,file_list));

%% Setup codes
if ~get_iscromer
    code_shift = 3;
    code_sample_on = 27;
    code_test = 29;
    code_short = 4;
    code_long = 30;
else
%     % Old method ... this one catches some bugs
%     code_shift = 3;
%     code_sample_on = 23;
%     code_test = 25;
%     code_short = 175;
%     code_long = 26;
    
    code_shift = 4;
    code_sample_on = 23;
    code_test = 25;
    code_short = 176;
    code_long = 27;
end

%% Add match nonmatch trial info to md

for i = 1:length(md)
    [md(i).sotr_match,md(i).sotr_nonmatch] = decode_match_nonmatch(md(i));
end

[aet,ae,sot,etr,sotr,match,nonmatch] = struct_pull(md,{'all_encode_times','all_encodes','sample_on_times','each_trial_response','sample_on_trials','sotr_match','sotr_nonmatch'});


%% Extract match nonmatch info and traces
plot_debug = 1;
if plot_debug; figure; end
for i = 1:length(ae)
    sample_on_encodes = find(ae{i} == code_sample_on);
    ind_err = etr{i} ~= 0;
    ind_err = ind_err(sotr{i});

    inds = sample_on_encodes(match{i}(sotr{i}) & ~ind_err);
    indl = sample_on_encodes(nonmatch{i}(sotr{i}) & ~ind_err);
%     ind_rest = ind27(ind_rest);
    ind_err = sample_on_encodes(ind_err);


    Npoints = 16;
    [indg] = get_centered_encodes(inds,Npoints);
    [aesq, tcent] = lookup_plotting_data(indg,inds,ae{i},aet{i});
    if plot_debug; clf; h1=plot(tcent,aesq,'b'); end
    ts{i} = tcent; aesegs{i} = aesq;

    [indg] = get_centered_encodes(indl,Npoints);
    [aesq, tcent] = lookup_plotting_data(indg,indl,ae{i},aet{i});
    if plot_debug; hold on; h2=plot(tcent,aesq,'r'); end
    tl{i} = tcent; aesegl{i} = aesq;
    
    [indg] = get_centered_encodes(ind_err,Npoints);
    [aesq, tcent] = lookup_plotting_data(indg,ind_err,ae{i},aet{i});
    if plot_debug; hold on; h2=plot(tcent,aesq-100,'g'); end
    terr{i} = tcent; aesegerr{i} = aesq;
    if plot_debug; 
        i
        pause(1);
    end
    
end

%% Plot all data

t = horzcat(ts{:});
x = horzcat(aesegs{:});
figure; plot(t,x,'b.');

t = horzcat(tl{:});
x = horzcat(aesegl{:});
hold on; plot(t,x,'r.');

t = horzcat(terr{:});
x = horzcat(aesegerr{:});
hold on; plot(t,x-100,'g.');


%% Plot anomalies



