



function [A,file_list,textdata] = load_conds(path, plot_on)
%% out = load_images_groups(path)

    if nargin < 2
        plot_on = 0;
    end
    
    file_list = get_filelist(path,'*.con');
    
    % Search through all files
    for i = 1:length(file_list)
        filename = fullfile(path,file_list{i});
        
        DELIMITER = ' ';

        % Import the file
        newData1 = importdata(filename, DELIMITER, 1);
        newData2 = importdata(filename, DELIMITER, 2);
        
        A{i} = [[newData1.data(1,:) NaN]; newData2.data(1:end,:)];
        textdata{i} = newData1.textdata;
        
    end
    
    if plot_on
        
        %% Plot all codes
        fprintf('Plotting some conds files... (hit any key)\n');
        pause;
        for i = 1:5;
            imagesc(A{i})
            fprintf('Hit any key \n')
            pause
        end
        
        %% Plot unique sample test1 and test2 codes...
        fprintf('Plotting unique sample test1 and test2 codes... (hit any key)\n');
        pause;
        usamples = calc_uniques(A,4);
        utest1 = calc_uniques(A,5);
        utest2 = calc_uniques(A,6);
%         figure;
%         subplot(221); imagesc(usamples); ylabel('File #'); xlabel('Array of unique codes'); title('Unique sample codes');
%         subplot(222); imagesc(utest1); ylabel('File #'); xlabel('Array of unique codes'); title('Unique test1 codes');
%         subplot(223); imagesc(utest2); ylabel('File #'); xlabel('Array of unique codes'); title('Unique test2 codes');
        
        figl;
        subplot(131); plot(unique(usamples)); ylabel('Unique codes in all files'); title('Sample Codes');
        subplot(132); plot(unique(utest1)); ylabel('Unique codes in all files'); title('Test 1 Codes');
        subplot(133); plot(unique(utest2)); ylabel('Unique codes in all files'); title('test 2 Codes');
        % Across all files, Sample codes are 201-234; test codes are <201
        % and >234.
        
    end
    
    %dlmread
end


function ucol2 = calc_uniques(A,col)

    ucol = cellfunu(@(x) unique(x(:,col)),A); lengths = cellfun(@length,ucol);
    
    ucol2 = zeros(length(ucol),max(lengths));
    for i = 1:length(ucol)
        ucol2(i,1:lengths(i)) = ucol{i};
    end
end

% 
% 
% function [A,file_list] = load_conds(path)
% %% out = load_images_groups(path)
% % Old version of function that reads through file line by line
% 
%     file_list = get_filelist(path,'*.con');
%     
%     % Search through all files
%     for i = 1:length(file_list)
%         filename = fullfile(path,file_list{i});
%         fid = fopen(filename);
%         
%         % Search through all lines in each file
%         j = 0;
%         tline = 'asdf';
%         while tline ~= -1
%             j=j+1;
%             tline = fgetl(fid);     % This is the command we want!
%             A{i}{j} = tline;
%         end
%         
%         fclose(fid);
%     end
%     
%     
% end