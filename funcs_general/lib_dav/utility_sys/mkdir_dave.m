

function mkdir_dave (output_path,suppress_output)
    % Deletes folder contents then recreates

%     fprintf(['Deleting all data in folder: \n' output_path '.\nPress any key to continue or CTRL-C to abort.\n']);
%     pause

    if nargin < 2
        suppress_output = 1;
    end

    if ~exist(output_path,'file')
        fprintf('Creating %s \n', output_path);
        %system( ['mkdir ' output_path]);
        [s1 s2] = mkdir(output_path);   % Outputs are supplied simply to suppress warning
                                        % message for existing folder.
    else
        if ~suppress_output; fprintf('Folder %s already exists. Doing nothing.\n',output_path); end
    end

    
end