

function bufferpath = get_buffer_home

    currbase = getpath('basepath');
    
    t1 = pwd;
    cd (currbase)
    currbase_resolved = pwd;        % Resolves symbolic links
    cd (t1);
    
    while (~strcmp(currbase_resolved,fileparts(pwd)) && ~strcmp(pwd,'/'))
        cd ..
    end
    
    [a b c] = fileparts(pwd);
    bufferpath = fullfile(getpath('path_buffer'),b);
    
    cd (t1);
    
end