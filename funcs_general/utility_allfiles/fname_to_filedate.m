

function ID = fname_to_filedate (fname)
%   generates a unique numeric identifier for each file name string
%   this unique Id likely contains the file date

    fname2 = fname(2:7);
    ID=str2num(fname2);
    
    if strcmp(fname(1),'L')
        ID=ID+800000;
    elseif strcmp(fname(1),'O')
        ID=ID+900000;
    else
        error('Monkey number not found');
    end
end