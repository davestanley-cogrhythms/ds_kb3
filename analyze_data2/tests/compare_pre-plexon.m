




datapath =('/Users/davestanley/Other/data/roy/CRC_Miller_data/intmatrix');
datapath_pre = ('/Users/davestanley/Other/data/roy_preplexon/CRC_Miller_data/intmatrix');
metadatapath = '/Users/davestanley/Other/data/roy/ctgX_merged_final/sansunits';


file_list = get_filelist(datapath_pre);
Nfiles = length(file_list);

% figure('Position',[  1          26        1280         769]);

file_range = 1:Nfiles;
%file_range = 2:Nfiles;

for i = file_range
    
    filename = file_list{i};
    
    d = load(fullfile(datapath,filename));
    dp = load(fullfile(datapath_pre,filename));
    md = load(fullfile(metadatapath,filename));
    
    elecnum = 1;
    lfp_name =md.lfp_names{elecnum}; lfp_name = strrep(lfp_name,'_',' ');
    x=double(d.lfpmat_int16(:,elecnum));
    xp=double(dp.lfpmat_int16(:,elecnum));
    
    sf = 100; %Shorten factor;
    
    xs=x(1:round(length(x)/sf),1);
    xps=xp(1:round(length(xp)/sf),1);
    
    fs = 1e3;
    t=1:length(xs);t=t/fs;
    
    figl;
    mylims = [50 70];
    subplot(2,2,1);
    plot(t,xps,'k'); legend('pre-plexon'); xlim(mylims)
    subplot(2,2,2); hold on; plot(t,xs,'b'); legend('post-plexon'); xlim(mylims)
    subplot(2,2,[3 4]); plot(t,xps,'k'); hold on; plot(t,xs,'b'); legend('pre-plexon','post-plexon'); xlim([57 60])

    subplot(2,2,1); title([filename ' ' lfp_name]);
    
end
