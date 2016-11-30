

function [sp_out] = map_sp(mode1, mode2,mypairs1,mypairs2,sp,md,bad_any1,bad_any2,sp_threshold)
    % sp is assumed to be in the form of mode1
    % we want to convert it to be of form of mode2
    
    if ~exist('sp_threshold','var')
        sp_threshold = 10;
    end
    
    plot_debug = 0;

    %data_type1 = mode2datatype(mode1);
    %data_type2 = mode2datatype(mode2);
    
    [pairs_type1, pairs_mode1 ]= mode2pairstype(mode1);
    [pairs_type2, pairs_mode2 ]= mode2pairstype(mode2);
    
    fprintf(['Converting ' pairs_type1 ' to ' pairs_type2 '\n']);
    
    if pairs_mode1 == pairs_mode2
        %% Same data types
        sp_out = sp;
    elseif ( pairs_mode1 == 0 || pairs_mode1 == 2 || pairs_mode1 == 7 || pairs_mode1 == 6) && pairs_mode2 == 4     % Units/PSD to FFC
        %%
        conversion_mode = 2;
%         tic
        switch conversion_mode
            case 1
                % Deleted this method...

            case 2

                % % % % Second method - more realistic % % % % %
                    % For each pair, look at all units associted with both
                    % electrodes. If any are true, mark as true (liberal). In other
                    % words, electrode pair is considered significant if any of its
                    % electrode pairs are significant.
                if pairs_mode1 ~= 6
                    ulpairs = build_unit_LFP_pairs_all(md);
                    sp_index_2_elect = ulpairs(2,:);
                else
                    sp_index_2_elect = mypairs1(:,1);
                end

                clear inds
                sp_out = false(size(mypairs2,1),size(sp,2));
                for i = 1:size(mypairs2,1)
                    p1 = mypairs2(i,1); p2=mypairs2(i,2);
                    inds{i} = find(sp_index_2_elect == p1 | sp_index_2_elect == p2);     % indices of units associated with the ith pair
                    sp_curr = sp(inds{i},:);
                    sp_out(i,:) = any(sp_curr,1);                                % Uncomment for either electrode carrying sig. unit.
                    %sp_out(i,:) = all(sp_curr,1) * double(~isempty(sp_curr));   % Uncomment for both electrodes carrying sig. Unit.
                                                    % ^^ Note - for some all returns trues if sp_curr is empty. This
                                                    %    forces the entries to be zero.
                end
        end
%         toc

    elseif ( pairs_mode1 == 0 || pairs_mode1 == 2 || pairs_mode1 == 7) && pairs_mode2 == 6     % Units to PSD
        %%
        ulpairs = build_unit_LFP_pairs_all(md);
        sp_index_2_elect = ulpairs(2,:);

        clear inds
        sp_out = false(size(mypairs2,1),size(sp,2));
        for i = 1:size(mypairs2,1)
            p1 = mypairs2(i,1);
            inds{i} = find(sp_index_2_elect == p1);     % indices of units associated with the ith pair
            sp_curr = sp(inds{i},:);
            sp_out(i,:) = any(sp_curr,1);
        end
    elseif pairs_mode1 == 6 && pairs_mode2 == 4     % PSD to FFC
        %%
        unit2elect = mypairs2(:,1);
        error('Incomplete');
    elseif pairs_mode1 == 4 && ( pairs_mode2 == 0 || pairs_mode2 == 2 || pairs_mode2 == 7 || pairs_mode2 == 6)     % FFC to Units/PSD
        %%
        
            % For each unit, determine which ensembles it's part of.
            % I.e. find all electrode pairs including it's electrode.
            % Assign identity accordingly
        
        % Parameters
        take_bads_into_account = 0;     % 0-Consider total true electrodes as a percentage of all electrodes; 1-As 0, but percentage of all GOOD electrodes
        if pairs_mode2 ~= 6
            ulpairs = build_unit_LFP_pairs_all(md);
            units2elect = ulpairs(2,:);
            mypairs2_elects = units2elect((mypairs2(:,1)));
        else
            mypairs2_elects = mypairs2(:,1);
        end
        
        
        %sp_out = false(size(mypairs2,1),size(sp,2));
        sp_trues = zeros(size(mypairs2,1),size(sp,2));
        pairs_per_electrode = zeros(size(mypairs2,1),1);
        for i = 1:length(mypairs2_elects)
            curr_elect = mypairs2_elects(i);
            ind = (mypairs1(:,1) ==  curr_elect | mypairs1(:,2) == curr_elect);
            sp_curr = sp(ind,:);
            %sp_out(i,:) = any(sp_curr,1);
            sp_trues(i,:) = sum(sp_curr,1);
            if ~take_bads_into_account
                pairs_per_electrode(i) = sum(ind);
            else
                pairs_per_electrode(i) = sum(~bad_any1(ind));   % Take into account "bad electrodes", which assume cannot possibly be true.
            end
        end
        
        percent_trues = sp_trues ./ repmat(pairs_per_electrode(:),[1,size(sp,2)])* 100;
        sp_out = percent_trues > sp_threshold; % Set this to 0 for the equivalent of "any" sp's being true.
                                                % This is essentially the fraction
                                                % of the unit's corresponding electrode
                                                % pairs that must have sp=1 to be
                                                % considered a significant unit.
        
        if plot_debug
            figure;
            plot(percent_trues)
        end
    elseif pairs_mode1 == 6 && ( pairs_mode2 == 0 || pairs_mode2 == 2 || pairs_mode2 == 7)      % PSD to Units
        %%
        ulpairs = build_unit_LFP_pairs_all(md);
        units2elect = ulpairs(2,:);
        mypairs2_elects = units2elect((mypairs2(:,1)));
        
        
        sp_out = false(size(mypairs2,1),size(sp,2));
        for i = 1:length(mypairs2_elects)
            curr_elect = mypairs2_elects(i);
            ind = find(mypairs1(:,1) ==  curr_elect);
            if length(ind) > 1; error('A single unit was found on more than one electrode. Something wrong.'); end
            sp_out(i,:) = sp(ind,:);
        end
        
    elseif pairs_mode1 == 4 && pairs_mode2 == 3                     % FFC to unitLFP pairs
        %%
        
            % For each unit, determine which ensembles it's part of.
            % I.e. find all electrode pairs including it's electrode.
            % Assign identity accordingly
        
        % Parameters
        take_bads_into_account = 0;     % 0-Consider total true electrodes as a percentage of all electrodes; 1-As 0, but percentage of all GOOD electrodes
        ulpairs = build_unit_LFP_pairs_all(md);
        units2elect = ulpairs(2,:);
        mypairs2_elects = [units2elect((mypairs2(:,1)))'  mypairs2(:,2)];
        
        %sp_out = false(size(mypairs2,1),size(sp,2));
        sp_trues = zeros(size(mypairs2,1),size(sp,2));
        pairs_per_electrode = zeros(size(mypairs2,1),1);
        for i = 1:length(mypairs2_elects)
            curr_elect = mypairs2_elects(i);
            ind = (mypairs1(:,1) == mypairs2_elects(i,1) & mypairs1(:,2) == mypairs2_elects(i,2)) | ...
                (mypairs1(:,2) == mypairs2_elects(i,1) & mypairs1(:,1) == mypairs2_elects(i,2));
            %warning('Check this works');
            sp_curr = sp(ind,:);
            %sp_out(i,:) = any(sp_curr,1);
            sp_trues(i,:) = sum(sp_curr,1);
            if ~take_bads_into_account
                pairs_per_electrode(i) = sum(ind);
            else
                pairs_per_electrode(i) = sum(~bad_any1(ind));   % Take into account "bad electrodes", which assume cannot possibly be true.
            end
        end
        
        sp_out = sp_trues > 0; 
       
        if plot_debug
            figure;
            plot(percent_trues)
        end
        
    elseif ( any(pairs_mode1 == [0,2]) && pairs_mode2 == 7 ) || ( pairs_mode1 == 7 && any(pairs_mode2 == [0,2]) ) % SFC to units or units to SFC
        sp_out = sp;
    else
        %%
        error('Unknown combination of pair types');  
    end
end

