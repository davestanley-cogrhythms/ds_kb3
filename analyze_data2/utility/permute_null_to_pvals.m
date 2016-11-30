



% sfc_mode = 22.40151111; permute_null_to_pvals(sfc_mode,3);
% sfc_mode = 22.40171111; permute_null_to_pvals(sfc_mode,2);
% sfc_mode = 22.45171111; permute_null_to_pvals(sfc_mode,3);
%sfc_mode = 22.44151110; permute_null_to_pvals(sfc_mode,2);
% sfc_mode = 22.451511113; permute_null_to_pvals(sfc_mode,2); permute_null_to_pvals(sfc_mode,3);
% sfc_mode = 2.20171110; permute_null_to_pvals(sfc_mode,2);
% sfc_mode = 2.00151110; permute_null_to_pvals(sfc_mode,3);
% sfc_mode = 41.6514111; permute_null_to_pvals(sfc_mode,3); permute_null_to_pvals(sfc_mode,2);
% sfc_mode = 41.6515111; permute_null_to_pvals(sfc_mode,3); permute_null_to_pvals(sfc_mode,2);
% sfc_mode = 41.6015111; permute_null_to_pvals(sfc_mode,3); permute_null_to_pvals(sfc_mode,2);


% sfc_mode = 52.70020100; permute_null_to_pvals(sfc_mode,4);
% sfc_mode = 52.74020100; permute_null_to_pvals(sfc_mode,4);
% sfc_mode = 52.74030100; permute_null_to_pvals(sfc_mode,4);
% sfc_mode = 52.75020100; permute_null_to_pvals(sfc_mode,4);
% sfc_mode = 52.75030100; permute_null_to_pvals(sfc_mode,4);

% sfc_mode = 52.700101; permute_null_to_pvals(sfc_mode,4);
% sfc_mode = 52.740101; permute_null_to_pvals(sfc_mode,4);
% sfc_mode = 52.750101; permute_null_to_pvals(sfc_mode,4);
% sfc_mode = 23.4014111; permute_null_to_pvals(sfc_mode,4);

% sfc_mode = 52.740201; permute_null_to_pvals(sfc_mode,4)
% sfc_mode = 52.740301; permute_null_to_pvals(sfc_mode,4)
% sfc_mode = 52.750201; permute_null_to_pvals(sfc_mode,4)
% sfc_mode = 52.750301; permute_null_to_pvals(sfc_mode,4)

function permute_null_to_pvals(sfc_mode,curr_stage)


    %% Setup vars
    if ~exist('sfc_mode','var'); sfc_mode=22.4415111; end
    if ~exist('curr_stage','var'); curr_stage=3; end

    %% First, get list of files to fix

    [~, mode_subgroups] = decode_sfc_mode(sfc_mode);
    [fname_suffix] = build_sfcmode(sfc_mode, mode_subgroups);
    outpath = fullfile(getpath('path_buffer_curr'),['test3mode_' num2str(sfc_mode,12) fname_suffix],['stage_' num2str(curr_stage)]);
    outpath2 = fullfile(outpath,'trimmed');
    mkdir_dave (outpath2);


    fn = dir(fullfile(outpath,'*.mat'));
    fn = {fn.name};


    %% Perform the operation
    parfor f = 1:length(fn)
        myconvert(outpath,fn,f,outpath2)
    end
    
    
    %% Rearrange things
%     keyboard
    system(['mkdir ' fullfile(outpath,'full')]);
    system(['mv ' fullfile(outpath,'*.mat') ' ' fullfile(outpath,'full')]);
    system(['mv ' fullfile(outpath2,'*.mat') ' ' fullfile(outpath)]);
    system(['rmdir ' fullfile(outpath2)]);
    

end


function [pvals, critval]= calcp(Cave1,Cave2,Cavep1,Cavep2, do_phi)

    % Get diff
    dCave = Cave1 - Cave2;
    if do_phi; dCave = get_min_phi_diff(dCave); end
    dCave = abs(dCave);

    % Git diff null
    dCave_null = Cavep1 - Cavep2;
    if do_phi; dCave_null = get_min_phi_diff(dCave_null); end
    sz = size(dCave_null);
    
    % Centre both values about zero (new code)
    null_means = mean(dCave_null,3);
    dCave = dCave - null_means;
    dCave_null = dCave_null - repmat(null_means,[1,1,sz(3)]);

    
    %pvals = sum(dCave_null > abs(repmat(dCave,[1,1,sz(3)])),3) ./ (sz(3)*ones([sz(1:2)])); % One tailed
    pvals = sum( (dCave_null > abs(repmat(dCave,[1,1,sz(3)]))) | (dCave_null < -abs(repmat(dCave,[1,1,sz(3)]))) ,3) ./ (sz(3)*ones([sz(1:2)])); % Two tailed
    critval = prctile(dCave_null,95,3);
            
end
    

function z = get_min_phi_diff(z)
    
    ind=z>pi; z(ind)=z(ind)-2*pi; ind=z<-pi; z(ind)=z(ind)+2*pi;        % Bounds of maximum shift should be -pi to pi
    %z = abs(z);
    
end

function [zsc] = calc_normalized_dCave(Cave1,Cave2,Cavep1,Cavep2, do_phi)
 
    % Get diff
    dCave = Cave1 - Cave2;
    if do_phi; dCave = get_min_phi_diff(dCave); end
    dCave = abs(dCave);

    % Git diff null
    dCave_null = Cavep1 - Cavep2;
    if do_phi; dCave_null = get_min_phi_diff(dCave_null); end
    sz = size(dCave_null);
    
    % Centre both values about zero (new code)
    null_means = mean(dCave_null,3);
    dCave = dCave - null_means;
    dCave_null = dCave_null - repmat(null_means,[1,1,sz(3)]);

    zsc = (dCave)./mean(abs(dCave_null),3);
    
    %pvals = sum(dCave_null > abs(repmat(dCave,[1,1,sz(3)])),3) ./ (sz(3)*ones([sz(1:2)])); % One tailed
    pvals = sum( (dCave_null > abs(repmat(dCave,[1,1,sz(3)]))) | (dCave_null < -abs(repmat(dCave,[1,1,sz(3)]))) ,3) ./ (sz(3)*ones([sz(1:2)])); % Two tailed
    critval = prctile(dCave_null,95,3);
            
end

function myconvert(outpath,fn,f,outpath2)
    clear sfc
    
    % Check if output file already exists
    outname = fullfile(outpath2,fn{f});
    if exist(outname,'file')
        fprintf(['Output file exists..skipping: ' outname '\n']);
        return
    end
    
    
    % Load file
    fprintf(['Loading file ' fullfile(outpath,fn{f}) '\n'])
    temp = load(fullfile(outpath,fn{f}));
    sfc = temp.sfc;

    for i = 1:length(sfc)

        s = sfc{i};
        if isfield(s,'Cave1');
            s.pvals_Cave = calcp(s.Cave1,s.Cave2,s.Cavep1,s.Cavep2,0);
            if isfield(s,'phi1');
                s.pvals_phi = calcp(s.phi1,s.phi2,s.phip1,s.phip2,1);
            end


            s.zsc_Cave = calc_normalized_dCave(s.Cave1,s.Cave2,s.Cavep1,s.Cavep2,0);
            if isfield(s,'phi1');
                s.zsc_phi = calc_normalized_dCave(s.phi1,s.phi2,s.phip1,s.phip2,1);
            end


            if isfield(s,'phi1');
                s = rmfield(s,{'Cavep1','Cavep2','phip1','phip2'});
            else
                s = rmfield(s,{'Cavep1','Cavep2'});
            end
        end

        sfc{i} = s;
    end
    save(outname,'sfc');
end