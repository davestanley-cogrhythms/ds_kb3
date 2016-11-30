function ctgsetli_maxes = calc_ctgsetli_maxes(ctgsetli_good,ctgsall,coh_debias_mode)
    %% function ctgsetli_maxes = calc_ctgsetli_maxes(ctgsetli_good,ctgsall,coh_debias_mode)
    % coh_debias_modee code from beginning of file:
    % sfc_mode = *.??????X -> Correct for uneven trial number bias
    %                         1-Min trials across pairs of i0s single cell; 2-Min trials across pairs of i0s all cells
    %                         3-Min trials across all i0s single cell; 4-Min trials all i0s all cells
    
    Nctgs = size(ctgsetli_good,2);
    sumctgs = sum(ctgsetli_good);
    
    % Make sure Nctgs is even; 
            if mod(Nctgs,2) ~= 0
                %warning('The paried coh_debias_modee mode requires an length(i0) to be even. Pairs of ctgs must be compared.');
                % ACtually, it's okay - if it's odd we'll just leave the
                % end one alone
            end

    if ~isempty(ctgsall)
        ctgsetli_mins_allcells = ctgsall.ctgsetli_mins_allcells;
    end
    
    switch coh_debias_mode
        case 0
            %ctgsetli_maxes = sum(ctgsetli_good);    % No change
            ctgsetli_maxes = zeros(1,size(ctgsetli_good,2));  % Set all to zero to code for no change
        case 1
            ctgsetli_maxes = sumctgs;
            Npairs = floor(Nctgs/2);
            for pi = 1:Npairs
                inds = [pi*2-1, pi*2];
                ctgsetli_maxes(inds) = min(sumctgs(inds));
            end
        case 2
            ctgsetli_maxes = ctgsetli_mins_allcells;
            Npairs = floor(Nctgs/2);
            for pi = 1:Npairs
                inds = [pi*2-1, pi*2];
                ctgsetli_maxes(inds) = min(ctgsetli_mins_allcells(inds));
            end
        case 3
            ctgsetli_maxes = zeros(1,Nctgs);
            ctgsetli_maxes(:) = min(sumctgs);
        case 4
            ctgsetli_maxes = zeros(1,Nctgs);
            ctgsetli_maxes(:) = min(ctgsetli_mins_allcells);
    end
    
    % Make sure that we don't prune ctgsetli_maxes to be shorter than Ntrials_thresh
    if coh_debias_mode == 2 || coh_debias_mode == 4
        Ntrials_thresh = get_Ntrials_thresh;
        ctgsetli_maxes( ctgsetli_maxes < Ntrials_thresh ) = Ntrials_thresh;
    end
end