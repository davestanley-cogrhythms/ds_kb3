

function group = get_grouped_cells_and_data(group,sp,pls,f,mypairs,bad_any,plotmode,freqband_stats,funames1D,abscissa2,is_spectrogram,timeband_stats)

if ~exist('abscissa2','var'); abscissa2 = []; end
if ~exist('is_spectrogram','var'); is_spectrogram = 0; end
if ~exist('timeband_stats','var'); timeband_stats = []; end


%% Calculate group cells
% Extract units sensitive to given schemes
for i = 1:length(group)
    [group(i).cells, group(i).numcells]= get_grouped_cells(group(i),[sp]);
end

% Sanity check to make sure all bad sp's have already been removed
if size(sp,1) == length(bad_any)
    sum_goods = sum(sp);
    sum_goods_and_not_bad = sum(sp & repmat(~bad_any(:),1,size(sp,2)));
    if any(sum_goods ~= sum_goods_and_not_bad);
        warning('Bads are not being properly removed from sp calculation!');
    end
end

% Remove bad files
[group] = remove_bads_cells(bad_any,group);                 % Remove bad cells from cell lists


%% Groups Loading Setup
% Generate SFC data for plotting and statistics
for i = 1:length(group)
    switch plotmode
        case {1,2,4,5,6,7,8,9}
            [group(i).data, group(i).funames] = get_grouped_data(group(i).ctgs,group(i).cells,pls,funames1D);
            [group(i).metadata.mypairs] = get_grouped_data(ones(size(group(i).cells,2),1),group(i).cells,mypairs',funames1D)';
            abscissa = f;
            
            if ~is_spectrogram
                [group(i).datastats, group(i).freqband_stats] = calc_pls_stats(f,group(i).data,freqband_stats,'do_mean_ctgs',0);
            else
                [group(i).datastats, group(i).freqband_stats, group(i).timeband_stats] = calc_pls3D_stats(f,abscissa2,group(i).data,freqband_stats,timeband_stats,'ctg_singleton',1);         % Set ctg_singleton to 0 for pls; to 1 for group.data. This is because dimensionality is different!
            end
            
            if plotmode == 6 || plotmode == 7
                warning('This still might not work properly');
                myfunc = @(x) (fcomplex2real(mean(x))); Nbs = 100;
                group(i).data = bootstrp_and_scale(Nbs,myfunc,group(i).data')';
                group(i).datastats = bootstrp_and_scale(Nbs,myfunc,group(i).datastats')';
                clear myfunc Nbs
                
                % Bootstrap override
%                 group(i).data = fcomplex2real(mean(group(i).data,2));
%                 group(i).datastats = fcomplex2real(mean(group(i).datastats,2));
            end
            
%             if do_fisher
%                 group = group.myfunc (group, 'data', @fisher_trans);
%                 group = group.myfunc (group, 'datastats', @fisher_trans);
%             end

        case 3
            [group(i).data, group(i).funames] = get_grouped_data(group(i).ctgs,group(i).cells,ns,funames1D);
            t = [1:size(ns,1)]*get_dt;
            t = t - 2;  % Centre to sample on time
            abscissa = t;
            timeband = get_stagesir(curr_stage_sfc)/1e3;
            [group(i).datastats, group(i).freqband_stats] = calc_pls_stats(t,group(i).data,timeband,'do_mean_ctgs',1);
    end
    group(i).xdata=abscissa;
    group(i).xdata2=abscissa2;
end