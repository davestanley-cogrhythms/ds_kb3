


function data = slice_data_closest(data,abscissa,dimension,val)
    % data = slice_data_closest(data,abscissa,dimension,val)
    
    % Find desired index (e.g. index at which f = 16 Hz)
    index = find_closest(abscissa,mean(val));
    
    % Bring chosen dimension to front
    data=permute(data,[dimension,1:dimension-1,dimension+1:ndims(data)]);
    
    % Slice data long that dimension at desired index
    data=data(index,:,:,:,:,:,:,:,:);
    
    % Restore original ordering
    data=permute(data,[2:dimension,1,dimension+1:ndims(data)]);
    
    % Squeeze data
    % data = squeeze(data);
    
    
end