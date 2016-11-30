

function gr = group_enumerate_stages(gr)
    % Builds the group to compare sample and delay stages

    gr_all = gr.all;
    gr_alt = group_swap.alt (gr.all);
    N1=length(gr_all);
    N2=length(gr_alt);
    
    gr_sd(1:N1*N2) = Grp_SD;       % Group class for comparing sample and delay
    k=1;
    for i = 1:N1
        for j = 1:N2
            gr_sd(k) = group_AND(gr_all(i),gr_alt(j));
            %gr_sd(k).inheritObj(group_AND(gr_all(i),gr_alt(j)));   % Alternate method to inheret object. More ugly..
            
            gr_sd(k).coords_alt = gr_alt(j).coords;
            gr_sd(k).right_alt = gr_alt(j).right;
            gr_sd(k).wrong_alt = gr_alt(j).wrong;
            gr_sd(k).diff_alt = gr_alt(j).diff;
            
            k=k+1;
        end
    end
    
    gr.sd = gr_sd;
    %%
    % Decider2's
    clear gc
    i=1; ind = horzcat(gr_sd.right) == 2 & horzcat(gr_sd.right_alt) == 2; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr2-Pr2';
    i=2; ind = horzcat(gr_sd.right) == 2 & horzcat(gr_sd.right_alt) == 0; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr2-Pr0';
    i=3; ind = horzcat(gr_sd.right) == 0 & horzcat(gr_sd.right_alt) == 2; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr0-Pr2';
    i=4; ind = horzcat(gr_sd.right) == 0 & horzcat(gr_sd.right_alt) == 0; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr0-Pr0';
    gr.sd_P2 = gc;
    
    
    %% All 9 possibilities
    % Decider1's
    clear gc
    i=1; ind = horzcat(gr_sd.right) == 1 & horzcat(gr_sd.right_alt) == 1; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr1-Pr1';
    i=2; ind = horzcat(gr_sd.right) == 1 & horzcat(gr_sd.right_alt) == 0; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr1-Pr0';
    i=3; ind = horzcat(gr_sd.right) == 0 & horzcat(gr_sd.right_alt) == 1; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr0-Pr1';
    i=4; ind = horzcat(gr_sd.right) == 0 & horzcat(gr_sd.right_alt) == 0; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr0-Pr0';
    gr.sd_P1 = gc;
    
    % Decider SP2
    clear gc
    i=1; ind = horzcat(gr_sd.right) == 2 & horzcat(gr_sd.right_alt) == 2; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr2-Pr2';
    i=2; ind = horzcat(gr_sd.right) == 2 & horzcat(gr_sd.right_alt) == 1; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr2-Pr1';
    i=3; ind = horzcat(gr_sd.right) == 2 & horzcat(gr_sd.right_alt) == 0; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr2-Pr0';
    gr.sd_sP2 = gc;
    
    % Decider SD2
    clear gc
    i=1; ind = horzcat(gr_sd.right) == 2 & horzcat(gr_sd.right_alt) == 2; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr2-Pr2';
    i=2; ind = horzcat(gr_sd.right) == 1 & horzcat(gr_sd.right_alt) == 2; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr1-Pr2';
    i=3; ind = horzcat(gr_sd.right) == 0 & horzcat(gr_sd.right_alt) == 2; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr0-Pr2';
    gr.sd_dP2 = gc;
    
    %%
    % Decider >1
    clear gc
    i=1; ind = horzcat(gr_sd.right) >= 1 & horzcat(gr_sd.right_alt) >= 1; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr-Pr';
    i=2; ind = horzcat(gr_sd.right) >= 1 & horzcat(gr_sd.right_alt) == 0; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr-NP';
    i=3; ind = horzcat(gr_sd.right) == 0 & horzcat(gr_sd.right_alt) >= 1; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='NP-Pr';
    i=4; ind = horzcat(gr_sd.right) == 0 & horzcat(gr_sd.right_alt) == 0; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='NP-NP';
    gr.sd_PNP = gc;
    
    %%
    % Decider all
    clear gc
    i=1; ind = horzcat(gr_sd.right) == 2 & horzcat(gr_sd.right_alt) == 2; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr2-Pr2';
    i=2; ind = horzcat(gr_sd.right) == 2 & horzcat(gr_sd.right_alt) == 1; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr2-Pr1';
    i=3; ind = horzcat(gr_sd.right) == 2 & horzcat(gr_sd.right_alt) == 0; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr2-Pr0';
    i=4; ind = horzcat(gr_sd.right) == 1 & horzcat(gr_sd.right_alt) == 2; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr1-Pr2';
    i=5; ind = horzcat(gr_sd.right) == 1 & horzcat(gr_sd.right_alt) == 1; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr1-Pr1';
    i=6; ind = horzcat(gr_sd.right) == 1 & horzcat(gr_sd.right_alt) == 0; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr1-Pr0';
    i=7; ind = horzcat(gr_sd.right) == 0 & horzcat(gr_sd.right_alt) == 2; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr0-Pr2';
    i=8; ind = horzcat(gr_sd.right) == 0 & horzcat(gr_sd.right_alt) == 1; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr0-Pr1';
    i=9; ind = horzcat(gr_sd.right) == 0 & horzcat(gr_sd.right_alt) == 0; gc(i) = group_merge(gr_sd(ind)); gc(i).legend='Pr0-Pr0';
    gr.sd_Pall = gc;
    

    %% SD AvB
    clear gc
    i=0;
    i=i+1; gc(i) = group_AND(gr_all(3),gr_alt(3)); gc(i).legend='PrA-PrA';
    i=i+1; gc(i) = group_AND(gr_all(2),gr_alt(2)); gc(i).legend='PrB-PrB';
    i=i+1; gc(i) = group_AND(gr_all(3),gr_alt(4)); gc(i).legend='PrA-Pr0';
    i=i+1; gc(i) = group_AND(gr_all(4),gr_alt(3)); gc(i).legend='Pr0-PrA';
    i=i+1; gc(i) = group_AND(gr_all(2),gr_alt(4)); gc(i).legend='PrB-Pr0';
    i=i+1; gc(i) = group_AND(gr_all(4),gr_alt(2)); gc(i).legend='Pr0-PrB';
    i=i+1; gc(i) = group_AND(gr_all(4),gr_alt(4)); gc(i).legend='Pr0-Pr0';
    gr.sd_SpecAll = gc;
    
    
end

