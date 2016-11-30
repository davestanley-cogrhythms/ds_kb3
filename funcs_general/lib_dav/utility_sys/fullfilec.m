
function out = fullfilec(varargin)
%     Like fullfile, but can accept a cell array as the last argument. If the last
%     argument is an ordinary string, then it will just return the same output as
%     fullfile.
%     USAGE
%     f = fullfile(folderName1, folderName2, ..., fileName)
%     f = fullfile(folderName1, folderName2, ..., cell_array_of_filenames)


    if iscellstr(varargin{end})
        out = cellfun(@(x) fullfile(varargin{1:end-1},x),varargin{end},'UniformOutput',0);
    else
        out = fullfile(varargin{:});
    end
    
end
