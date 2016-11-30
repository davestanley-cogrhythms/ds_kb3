function userdir = getpath_home

    if ispc; userdir= getenv('USERPROFILE'); 
    else; userdir= getenv('HOME'); 
    end

end