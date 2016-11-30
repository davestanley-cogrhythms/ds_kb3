
function out = zeros_gpuwrapper(varargin)
    % Usage = use like normal zeros command, except supply some sample data
    % at the end. If this sample data is of type gpuArray, this command
    % will generate zeros on the GPU instead of CPU.
    
    argin = varargin(1:end-1);
    sample_data = varargin{end};
    
    if is_gpuArray(sample_data)
        out = zeros(argin{:},'gpuArray');
    else
        out = zeros(argin{:});
    end
    
end