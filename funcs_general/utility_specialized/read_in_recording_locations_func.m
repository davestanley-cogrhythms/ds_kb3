


function [el_struct] = read_in_recording_locations_func(do_monkeyL)

plot_on_PFC = 0;
plot_on_ITC = 0;
plot_on_ITC_3d = 0;


rotate_flag = 0;

PathName = getpath('path_electrode_locations');
if do_monkeyL
    FileName = 'LUCY right well recordings catdog.xls';
else
    FileName = 'Ozzy right PFC recording sites.xls';
end

[~,~,whole_thing] = xlsread(fullfile(PathName,FileName));

if ~isempty(strfind(FileName,'LUCY'))
    whole_thing = whole_thing(6:end,:);
end

if ~isempty(strfind(FileName,'Ozzy'))
    whole_thing = whole_thing(3:end,:);
end

% Pull out date information
dateinfo = cell2mat(whole_thing(:,2));
electrode_num = cell2mat(whole_thing(:,3));

% Search & replace nans in dateinfo
dateinfo_old = dateinfo;
for i = 1:length(dateinfo)
    if ~isnan(dateinfo(i))
        currdate = dateinfo(i);
    else
        dateinfo(i) = currdate;
    end
end

%%
%  divide array up  PFC
% convert from cell to str to integer
letters = cell2mat(whole_thing(:,4));    
if ischar(letters)
    letters = char(letters);
    letters = double(letters) - 97;   
end



AP = letters(:,1) ;             % This is kind of redundant
ML = cell2mat(whole_thing(:,5)) ;

% Add some noise
% noise_mag = .5;
% AP = AP + unifrnd(-noise_mag,noise_mag,size(AP));
% ML = ML + unifrnd(-noise_mag,noise_mag,size(AP));

%
unit_or_not = cell2mat((whole_thing(:,7)));
unit_or_not(isnan(unit_or_not)) = 0;

unit_locations = find(0 < unit_or_not & unit_or_not < 999);
bone_locations = find(unit_or_not == 999);

rotate_by = cell2mat(whole_thing(:,6));
rotate_by(isnan(rotate_by)) = 0;

if rotate_flag == 1
new_AP = AP;
new_ML = ML;

new_AP(rotate_by~=0) =  AP(rotate_by~=0) .* cos((rotate_by(rotate_by~=0)*pi/180)) - ML(rotate_by~=0) .* sin((rotate_by(rotate_by~=0)*pi/180));
new_ML(rotate_by~=0) =  ML(rotate_by~=0) .* cos((rotate_by(rotate_by~=0)*pi/180)) + AP(rotate_by~=0) .* sin((rotate_by(rotate_by~=0)*pi/180));

else
    new_AP = AP;
new_ML = ML;
end

%  divide array up  ITC
% convert from cell to str to integer

ITC_AP = cell2mat(whole_thing(:,11));
ITC_AP(isnan(ITC_AP)) = [];

ITC_ML = cell2mat(whole_thing(:,12));
ITC_ML(isnan(ITC_ML)) = [];

ITC_dp = cell2mat(whole_thing(:,13));
ITC_dp(isnan(ITC_dp)) = [];

%%
% plot coordinates
if plot_on_PFC
    figure('Name', 'Recording locations PFC', 'Position', [300 300 600 600],'NumberTitle','off');

    title([FileName(1:end-4), ' PFC  as of ', date])

    hold on
    rectangle('Position',[0,0,25,25],'Curvature',[1,1])
    rectangle('Position',[3,3,19,19],'Curvature',[1,1], 'LineStyle','--')
    plot(new_AP, new_ML, 'ok', new_AP(unit_locations), new_ML(unit_locations), 'og', new_AP(bone_locations), new_ML(bone_locations), 'or')
    text(25, 12.5, '>', 'FontSize', 30)
    hold off

    legend('no units', 'units', 'bone')
    xlabel('Post <-------------------> Ant')
    ylabel('Lat <-------------------> Med')

    axis square
    axh = gca;
    set(axh,'XTick',[1:1:25]);
    set(axh,'YTick',[1:1:25]);

    grid
end

%%
if plot_on_ITC
    figure('Name', 'Recording locations ITC', 'Position', [200 200 600 600],'NumberTitle','off');

    title([FileName(1:end-4), ' ITC  as of ', date])

    hold on
    rectangle('Position',[-12.5,-12.5,25,25],'Curvature',[1,1])
    rectangle('Position',[-9.5,-9.5,19,19],'Curvature',[1,1], 'LineStyle','--')
    plot(ITC_AP, ITC_ML, 'ok')
    text(12.5, 0, '>', 'FontSize', 30)
    hold off

    legend('units')
    xlabel('Post <---------well zero = 39.3 ----------> Ant')
    ylabel('Lat <-------------------> Med')

    axis square
    %grid
end

%%
if plot_on_ITC_3d
    figure('Name', 'Recording locations ITC', 'Position', [200 200 600 600],'NumberTitle','off');

    title([FileName(1:end-4), ' ITC  as of ', date])

    hold on
    rectangle('Position',[-12.5,-12.5,25,25],'Curvature',[1,1])
    rectangle('Position',[-9.5,-9.5,19,19],'Curvature',[1,1], 'LineStyle','--')
    plot3(ITC_AP, ITC_ML, ITC_dp, 'ok')
    text(12.5, 0, '>', 'FontSize', 30)
    hold off

    line([-10  10],[0 0],[24.6 24.6], 'color', 'red')
    text([-10  10],[0 0],[25.3 25.3], 'TPO')
    line([-10  10],[0 0],[25.6 25.6], 'color', 'red')

    line([-10  10],[0 0],[26.9 26.9], 'color', 'red')
    text([-10  10],[0 0],[27 27], 'TEa')
    line([-10  10],[0 0],[27.04 27.04], 'color', 'red')

    line([-10  10],[0 0],[28.24 28.24], 'color', 'red')
    text([-10  10],[0 0],[28.5 28.5], 'TE1')

    legend('units')
    xlabel('Post / Ant')
    ylabel('Lat / Med')
    zlabel('depth')

    set(gca, 'ZDir', 'reverse')
    axis square
    view(-15, 35)
    %grid
end

%% Pack up data
el_struct.AP = new_AP;
el_struct.ML = new_ML;
el_struct.unit_locations = unit_locations;
el_struct.bone_locations = bone_locations;
el_struct.dateinfo = dateinfo;
el_struct.electrode_num = electrode_num;
  

end
