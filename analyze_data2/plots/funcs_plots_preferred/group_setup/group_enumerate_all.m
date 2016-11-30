

function [gr_SD_QP, gr_S_QP] = group_enumerate_all(gr_all,gr_single)
    % Builds the group to compare sample and delay stages

%     gr_sd = gr.sd;
%     gr_all = gr.all;
%     gr_alt = group_swap.alt (gr.all);
%     
%     gr_single = gr.single;
    
    N_ctgs_base = get_Nctgs_base;
    N_ctgs_extra = get_Nctgs_extra;
    N_pairs_extra = round(N_ctgs_extra / 2);
    
    % Sample deciders (P - preferred, paradigm)
    clear gc
    i=1; ind = horzcat(gr_all.right) >= 1; gc(i) = group_merge(gr_all(ind)); gc(i).legend='SP+';
    i=2; ind = horzcat(gr_all.right) == 0; gc(i) = group_merge(gr_all(ind)); gc(i).legend='SP0';
    gr_SP = gc;
    
    % Delay deciders
    clear gc
    gc = group_swap.alt(gr_SP);
    gc(1).legend = 'DP+';
    gc(2).legend = 'DP0';
    gr_DP = gc;
    
    % Sample cue deciders (Q)
    clear gc
    gc = gr_single;
    swap_map = [1 2; [9 10] ];             % schA schB
    swap_map_pref = unique(round(swap_map/2)','rows')';
    gc = group_swap.crit (gc,swap_map_pref); gc = group_swap.ctgs (gc,swap_map);
    gc(1).legend = 'SQ+';
    gc(2).legend = 'SQ0';
    gr_SQ = gc;
    
    % Delay cue deciders (Q)
    clear gc;
    gc = group_swap.alt(gr_SQ);
    gc(1).legend = 'DQ+';
    gc(2).legend = 'DQ0';
    gr_DQ = gc;
    
    
    % All combinations
    clear gc
    i=1; gc(i) = group_AND(gr_SP(1),gr_DP(1),gr_SQ(1),gr_DQ(1));
    i=2; gc(i) = group_AND(gr_SP(1),gr_DP(1),gr_SQ(1),gr_DQ(2));
    i=3; gc(i) = group_AND(gr_SP(1),gr_DP(1),gr_SQ(2),gr_DQ(1));
    i=4; gc(i) = group_AND(gr_SP(1),gr_DP(1),gr_SQ(2),gr_DQ(2));
    
    i=5; gc(i) = group_AND(gr_SP(1),gr_DP(2),gr_SQ(1),gr_DQ(1));
    i=6; gc(i) = group_AND(gr_SP(1),gr_DP(2),gr_SQ(1),gr_DQ(2));
    i=7; gc(i) = group_AND(gr_SP(1),gr_DP(2),gr_SQ(2),gr_DQ(1));
    i=8; gc(i) = group_AND(gr_SP(1),gr_DP(2),gr_SQ(2),gr_DQ(2));
    
    i=9; gc(i) = group_AND(gr_SP(2),gr_DP(1),gr_SQ(1),gr_DQ(1));
    i=10; gc(i) = group_AND(gr_SP(2),gr_DP(1),gr_SQ(1),gr_DQ(2));
    i=11; gc(i) = group_AND(gr_SP(2),gr_DP(1),gr_SQ(2),gr_DQ(1));
    i=12; gc(i) = group_AND(gr_SP(2),gr_DP(1),gr_SQ(2),gr_DQ(2));
    
    i=13; gc(i) = group_AND(gr_SP(2),gr_DP(2),gr_SQ(1),gr_DQ(1));
    i=14; gc(i) = group_AND(gr_SP(2),gr_DP(2),gr_SQ(1),gr_DQ(2));
    i=15; gc(i) = group_AND(gr_SP(2),gr_DP(2),gr_SQ(2),gr_DQ(1));
    i=16; gc(i) = group_AND(gr_SP(2),gr_DP(2),gr_SQ(2),gr_DQ(2));
    
    gr_SD_QP=gc;
    
    
    % Current stage; Cue + P
    clear gc
    i=1; gc(i) = group_AND(gr_SP(1),gr_SQ(1));
    i=2; gc(i) = group_AND(gr_SP(1),gr_SQ(2));
    i=3; gc(i) = group_AND(gr_SP(2),gr_SQ(1));
    i=4; gc(i) = group_AND(gr_SP(2),gr_SQ(2));
    gr_S_QP = gc;
    
end

