
function [sp,sr,ss,pr_fnames,pr_funames1D] = extract_spsr(stage, bad_any,use_bh_stepup)

    alpha = 0.05;
    
    if ~exist('use_bh_stepup','var');
        use_bh_stepup = 0; % Use BH Step Up algorithm to correct for multiple comparisons
    end
    
    if use_bh_stepup == 0
        warning('BH steputp is off!');
    end

    path_prefunits = fullfile(getpath('path_buffer_mypreferred'));
    pref_prefix = 'preferred_';    
    prefname = [pref_prefix 'stage_' num2str(stage)];   

    pref_struct = load(fullfile(path_prefunits,prefname));
    
    sp = vertcat(pref_struct.pschemes{:});
    sr = vertcat(pref_struct.spikerates{:});
    ss = vertcat(pref_struct.spikestd{:});
    pr_fnames = pref_struct.files_examined;
    pr_funames1D = [pref_struct.funames{:}];
    
    if ~use_bh_stepup
        sp = sp < 0.05/4;
    else
        %%
        cells1 = false(size(sp));
        cells2 = false(size(sp));
        cells3 = false(size(sp));
        clear critp1 critp2 cirtp3
        
%         % Control false discovery rate for all cells
%         for i=1:size(sp,2)
%             [cells1(:,i), critp1(i)]= bh_step1D(sp(:,i),alpha);
%         end
        
%         % Control false discovery rate for ONLY good cells. Since we're rejecting bad cells anyways, it makes sense to only consider good cells here.
%         tic
%         for i=1:size(sp,2)
%             [~, critp2(i)]= bh_step1D(sp(~bad_any,i),alpha);      % Just uses bad_any to get the crit value; there could still be bad cells in cells2. We remove these later.
%             cells2(:,i) = sp(:,i) <= critp2(i);
%         end
%         toc
        
        % bh stepup code 2D.
        
        if length(bad_any) ~= size(sp,1)
            warning('Size mismatch between number of cells and size bad_any. Setting all cells to good.');
            bad_any = false(1,size(sp,1));
        end
        
        
        tic
        [~, critp3]= bh_step(sp(~bad_any,:),alpha);
        cells3 = sp <= repmat(critp3,size(sp,1),1);
        toc
        %warning('Try replacing this code with the 2D bh_step function!');

        
%         find out why cells1 and cells 2 are different!!
%         [sum(cells1) ;sum(cells2); sum(cells3)]
        
        
        sp = cells3;
        %critp3
    end

end


