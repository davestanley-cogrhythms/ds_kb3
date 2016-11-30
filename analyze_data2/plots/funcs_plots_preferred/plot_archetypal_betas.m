
function plot_archetypal_betas(group,sfc,spc)
    % This code is used to find "archetypal" cells - i.e. those that fall
    % into one of the 16 groups and that --also-- have high beta SFC. The
    % high SFC cells are selected manually as follows:
%         1) Run the code with do_plott_ani0 = 1. This will cause the code
%         to cycle through all of the cells in each group, and you can
%         select which ones you want to use (generally those having high
%         beta SFC)
%         2) Remember the values you selected and plug them into the gsi
%         input in the test function. Now, when you run test group, it will
%         just plot that cell's SFC and will return you the file name and
%         information for that cell
    do_plott_ani0 = 0;
    

    % AdA, BdB group
    gr = group(1);
    chosen_cell = test_group (gr,sfc,spc,3,do_plott_ani0)
    
    
    % AdA, Bd0 group
    gr = group(9);
    chosen_cell = test_group (gr,sfc,spc,2,do_plott_ani0)
    chosen_cell = test_group (gr,sfc,spc,18,do_plott_ani0)
    
    % AdA, BdA group
%     gr = group(13);
%     chosen_cell = test_group (gr,sfc,spc,-1,do_plott_ani0)   % Can't find any cells
    
    
    % BdB, Ad0 group
    gr = group(3);
    chosen_cell = test_group (gr,sfc,spc,3,do_plott_ani0)
    chosen_cell = test_group (gr,sfc,spc,4,do_plott_ani0)
    
    % BdB, AdB group
%     gr = group(4);
%     chosen_cell = test_group (gr,sfc,spc,-1,do_plott_ani0)   % Can't find anything!
    
end

function chosen_cell = test_group (gr,sfc,spc,gsi,do_plott_ani)

    gr.coords       % Confirm we are using the correct coordinates
    
    cells = find(any(gr.cells,2)); %Cells in this group
    chosen_cell = cells(gsi);      %gsi is the index of the cells we have manually chosen
    
    if do_plott_ani; figure; plott_ani(squeeze(sfc(:,cells,:))); else
        figure; plot(squeeze(sfc(:,chosen_cell,:))); end
    
    
    [fnum, unum, curr_funame] = autorun_unum2fnum(chosen_cell)
    spc(chosen_cell,:)
    
    
end

