

plot_units_comparison = 0;
curr_stage=3;
addpath(genpath('./funcs_supporting_local'));

% Load preferred units
path_prefunits = fullfile(getpath('path_buffer_mypreferred'));
load(fullfile(path_prefunits,'preferred_units_dave'));

% Load SFC survey & format survey data
path_SFC = fullfile(getpath('path_bufffer_SFCsurvey'));
%load(fullfile(path_SFC,'betaSFC_survey2'));
load(fullfile(path_SFC,'betaSFC_survey1_partial'));

so = cat(2,survey_out{:});
so = strfind(so,'b9');  % Pull out all data labelled as 'b9'
so = ~cellfun(@isempty,so);
Nsurvey = length(so);
ui_SFCb = zeros(1,size(sch_raw,1));
ui_SFCb(1:Nsurvey) = so;


% Generate some indices for categorizing data
ui_anysig = or_across_rows(sch_raw); ui_anysig = (ui_anysig(:))';
ui_alldata = true(size(ui_anysig));
ui_bothtrue = (ui_SFCb == 1) & (ui_anysig == 1);
ui_outliers = (ui_SFCb == 1) & (ui_anysig == 0);

fprintf('%f percent of SFCb units are preferred units \n', 100*sum(ui_bothtrue)/sum(ui_SFCb));
fprintf('%f percent of SFCb units are not preferred units \n', 100*sum(ui_outliers)/sum(ui_SFCb));

% figure;
% hold on; plot(ui_anysig);
% hold on; plot(ui_SFCb,'kx')
% hold on; plot(ui(ui_bothtrue),ui_bothtrue(ui_bothtrue),'gx');
% hold on; plot(ui(ui_outliers),ui_outliers(ui_outliers),'rx');



file_range = 1:79;

% Load spiking traces
ns_struct = loadall_ra(5,5.1,'numspikes',[],file_range);
[ns, ns_fname, ns_uname] = merge_allfiles_ra(ns_struct); clear ns_struct


% Load SFC traces
sfc_struct = loadall_ra(curr_stage,2,'Cave',[],file_range);
[sfc, sfc_fname, sfc_uname] = merge_allfiles_ra(sfc_struct); clear sfc_struct

% Check that the correct files were loaded in both cases
if sum(~strcmp(sfc_uname,ns_uname)) ~= 0; fprintf(['File mismatch! Likely some files are missing for either '...
    'numspikes or sfc. Aborting \n']); end

if plot_units_comparison    % Don't plot this unless you want
    figure; plot_units_and_sfc(ns,sfc,ui_bothtrue,1)
    figure; plot_units_and_sfc(ns,sfc,ui_outliers,1)
    figure; plot_units_and_sfc(ns,sfc,ui_anysig,1)
    figure; plot_units_and_sfc(ns,sfc,uiSFCb,1)
    figure; plot_units_and_sfc(ns,sfc,ui_alldata,1)
end





