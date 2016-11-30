
function [in_ctg in_ctgnr] =  ctg_query(trial_number,md)
%     Input a trial number and metadata structure
%     Returns the relevant categories and non-relevant categories
%     associated with that trial number
%     INPUT
%     trial_number - 1 to # of trials
%     md - metadata array
%      
%     OUTPUT
%     in_ctg - 1x4 array with binary entries corresponding to ctg1-4. 1=is member; 0=not member
%     in_ctgnr - likewise for ctg_nr
    
    in_ctg = zeros(1,4);
    in_ctgnr = zeros(1,4);
    
    in_ctg(1) = sum(trial_number == md.ctg1_trials);
    in_ctg(2) = sum(trial_number == md.ctg2_trials);
    in_ctg(3) = sum(trial_number == md.ctg3_trials);
    in_ctg(4) = sum(trial_number == md.ctg4_trials);
    
    in_ctgnr(1) = sum(trial_number == md.ctg1_nr_trials);
    in_ctgnr(2) = sum(trial_number == md.ctg2_nr_trials);
    in_ctgnr(3) = sum(trial_number == md.ctg3_nr_trials);
    in_ctgnr(4) = sum(trial_number == md.ctg4_nr_trials);

%     ctg1_sampli = get_good_samples(md.ctg1_trials,md);
%     ctg2_sampli = get_good_samples(md.ctg2_trials,md);
%     ctg3_sampli = get_good_samples(md.ctg3_trials,md);
%     ctg4_sampli = get_good_samples(md.ctg4_trials,md);
%     
%     ctg1nr_sampli = get_good_samples(md.ctg1_nr_trials,md);
%     ctg2nr_sampli = get_good_samples(md.ctg2_nr_trials,md);
%     ctg3nr_sampli = get_good_samples(md.ctg3_nr_trials,md);
%     ctg4nr_sampli = get_good_samples(md.ctg4_nr_trials,md);
%     
    
    
end