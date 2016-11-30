

function conduct_sfc_survey(file_range)


sfc_mode = 2;
curr_stage = 2;
%datapath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode)]);
datapath = fullfile(getpath('path_buffer_curr'),['mode_' num2str(sfc_mode)],['stage_' num2str(curr_stage)]);
%datapath = fullfile(outpath,filename);

file_list = get_filelist(datapath);
Nfiles = length(file_list);

%figure('Position',[  1          26        1280         769]);
figure;
if  ~exist('file_range','var')
    file_range = 1:Nfiles;
    %file_range = [43];
end

for i = file_range
    filename = file_list{i};
    [survey_out{i}, units_examined{i}] = plot_cohfile(filename,datapath); 
    if strcmp(survey_out{i}{end},'q') || strcmp(survey_out{i}{end},'Q')
        break;
    end
    
end


files_examined=file_list;
fprintf('Press any key to save survey output');
pause
survey_path = fullfile(getpath('path_bufffer_SFCsurvey'),['stage' num2str(curr_stage)]);
mkdir_dave(survey_path)
save(fullfile(save_path,'betaSFC_survey2.mat'),'survey_out','files_examined','units_examined');

end


function [survey_out,unit_names_out] = plot_cohfile(filename,datapath) 
    take_survey = 1;
        subplot_on = 0;
    load (fullfile(datapath,filename))
    Nunits = size(sfc{9}.Cave,2);
    
    
    ctgs = [1 2 3 4 9];
    colourarr='bgrmk';
    
    if ~take_survey
        clf;
    end
    
    for j = 1:Nunits
        unit_names = sfc{1}.unit_names(j);
        unit_names=  strrep(unit_names,'_',' ');
        if ~take_survey
            subplotsq(Nunits,j);
            for i = 1:length(ctgs)
                hold on; plot(sfc{1}.f,sfc{ctgs(i)}.Cave(:,j),colourarr(i));
                title([filename ' ' unit_names])
                xlim([1 500])
                xlabel('freq (Hz'); ylabel('SFC');
            end
        else
            clf;
            for i = 1:length(ctgs)
                if subplot_on
                    subplotsq(length(ctgs),i);
                end

                hold on; plot(sfc{1}.f,sfc{ctgs(i)}.Cave(:,j),colourarr(i));
                title([filename ' ' unit_names])
                xlim([1 500])
                ylim([0 0.5])
                xlabel('freq (Hz)'); ylabel('SFC');
            end
        end
        
        if take_survey
            prompt = 'Type a comment or hit enter to continue q to quit: ';
            survey_out{j} = input(prompt,'s');
            
            if strcmp(survey_out{j},'q') || strcmp(survey_out{j},'Q')
                break;
            end
        else
            survey_out{j} = '';
        end
        unit_names_out{j} = unit_names;
        
    end
    
    if ~take_survey
        pause
    end

end

