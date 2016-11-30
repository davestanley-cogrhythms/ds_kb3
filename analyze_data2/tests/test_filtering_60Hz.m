


reload = 0;

datapath =('/Users/davestanley/Other/data/roy/CRC_Miller_data/intmatrix');
datapath_pre = ('/Users/davestanley/Other/data/roy_preplexon/CRC_Miller_data/intmatrix');
metadatapath = '/Users/davestanley/Other/data/roy/ctgX_merged_final/sansunits';


file_list = get_filelist(datapath_pre);
Nfiles = length(file_list);

% figure('Position',[  1          26        1280         769]);

file_range = 1:Nfiles;


if reload
    for i = file_range

        filename = file_list{i};

        %dc{i} = load(fullfile(datapath,filename));
        %dp{i} = load(fullfile(datapath_pre,filename));
        %mdc{i} = load(fullfile(metadatapath,filename));
    end
end


for i = file_range
    
    filename = file_list{i};
    
    d = dc{i};
    md = mdc{i};
    
    
    elecnum = 1;
    lfp_name =md.lfp_names{elecnum}; lfp_name = strrep(lfp_name,'_',' ');
    x=double(d.lfpmat_int16(:,elecnum));
    
    sf = 100; %Shorten factor;
    
    xs=x(1:round(length(x)/sf),1);

    fs = 1e3;
    t=1:length(xs);t=t/fs;
    
    xf = qif(t,xs,qif_notches([60 120],2));
    X=[xs xf];
    
    figl;
    %mylims = [50 70];
    mylims = [0 t(end)];
    hold on;
    clear fname;
    fname{1} = @plot;
    fname{2} = @plott_psd;
    %fname{3} = @plott_spect;
    plott_handles(fname,@subplotsq,t,X);
    %xlim(mylims)

    title([filename ' ' lfp_name]);
    
end
