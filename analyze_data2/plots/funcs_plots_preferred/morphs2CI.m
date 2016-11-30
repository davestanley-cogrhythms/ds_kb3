
function pls_CI = morphs2CI(ps,bad_any)
    %% pls_CI = morphs2CI(ps,bad_any)
    plot_debug = 0;
    
    
    calc_mode = 3;  % 1-My way of doing it - average 100,80,and 60% morphs to get WCD
                    % 2-Jefferson's way - only one-step differences to
                    %   get WCD/BCD
                    % 3-Jefferson's way - use one- and two-step differences
    
    
    do_abs_wcd = 1;     % 0-Only look at slope of WCD (not so good - this
                                % taking the average cancels out)
                        % 1-measures "roughness" of WCD  (Jefferson's)
                        
    

    switch calc_mode
        case 1
            % My way of doign it; average everything
            BCD = mean(ps(:,:,1:3),3) - mean(ps(:,:,5:7),3);
            BCD = abs(BCD);
            
%             BCD = cat(3,ps(:,:,1) - ps(:,:,7), ps(:,:,2) - ps(:,:,6),ps(:,:,3) - ps(:,:,5));
%             BCD = mean(abs(BCD),3);
            
            WCD = cat(3,ps(:,:,1) - ps(:,:,2), ps(:,:,2) - ps(:,:,3),ps(:,:,5) - ps(:,:,6),ps(:,:,6) - ps(:,:,7));
            if do_abs_wcd; WCD = mean(abs(WCD),3);
            else WCD = abs(mean(WCD,3)); end         % Note - does not requre sort_on = 1;
                                                     % Note - this is just equivalent to 100% morphs
                                                     %        minus 0% morphs
        %     % %
        case 2
            % One step differences
            BCD = abs(ps(:,:,3) - ps(:,:,5));
            WCD = cat(3,ps(:,:,1) - ps(:,:,2), ps(:,:,2) - ps(:,:,3),ps(:,:,5) - ps(:,:,6),ps(:,:,6) - ps(:,:,7));
            if do_abs_wcd; WCD = mean(abs(WCD),3);
            else WCD = abs(mean(WCD,3)); end         % Note - does not requre sort_on = 1;
                                                     % Note - this is just equivalent to 100% morphs
                                                     %        minus 0% morphs
            
            %WCD = (mean(WCD,3));
            
            
        case 3
            
            average_separately = 1;     % Should always set to 1.
            
            % One step differences
            BCD1 = abs(ps(:,:,3) - ps(:,:,5));
            WCD1 = cat(3,ps(:,:,1) - ps(:,:,2), ps(:,:,2) - ps(:,:,3),ps(:,:,5) - ps(:,:,6),ps(:,:,6) - ps(:,:,7));
             
            % Two step differences
            BCD2 = cat(3,ps(:,:,2) - ps(:,:,5), ps(:,:,3) - ps(:,:,6));
            WCD2 = cat(3,ps(:,:,1) - ps(:,:,3), ps(:,:,5) - ps(:,:,7));
            
            
            if average_separately
                % This way makes more sense - take average of one and two-step BCD and WCD results
                % rather than just lumping them together
                BCD1 = mean(abs(BCD1),3);
                WCD1 = mean(abs(WCD1),3);
                
                BCD2 = mean(abs(BCD2),3);
                WCD2 = mean(abs(WCD2),3);
                
                BCD = cat(3,BCD1,BCD2);
                WCD = cat(3,WCD1,WCD2);
                BCD = mean(abs(BCD),3);
                WCD = mean(abs(WCD),3);
            else
                BCD = cat(3,BCD1,BCD2);
                WCD = cat(3,WCD1,WCD2);
                BCD = mean(abs(BCD),3);
                WCD = mean(abs(WCD),3);
            end
            
            
            
            
    end
   
    
    
    pls_CI = (BCD - WCD) ./ (BCD + WCD);
    
    
    if plot_debug       % Plots bar graph (like Fig 9b roy)
        if ~exist('bad_any','var'); bad_any = false(size(ps,2),1); end
        figure; bar(squeeze(mean(ps(20,~bad_any,:),2)))
    end
    
    if 0
        % Testing code
        %%
        ind = 20;
        ps2 = squeeze(ps(ind,~bad_any,:));
        %figure; bar(mean(ps2,1))



        pls_CI2 = pls_CI(ind,~bad_any);

        [~,I] = sort(pls_CI2,'descend');

        %
        figure; plot(mean(ps2(I(1:100),:),1)); ylim([0.34 0.5]);
        yl = ylim;
        hold on; bar(mean(ps2(I(1:100),:),1)); ylim([0.34 0.5]);
        ylim (yl);

        figure; plot(mean(ps2(I(end-100:end),:),1)); ylim([0.34 0.5]);
        yl = ylim;
        hold on; bar(mean(ps2(I(end-100:end),:),1)); ylim([0.34 0.5]);
        ylim (yl);


        %%
        figure
        for i = 3000:10:length(I)
            clf;
            plot(ps2(I(i),:)); yl = ylim;
            hold on; bar(ps2(I(i),:));
            title(pls_CI2(I(i)));
            ylim (yl);
             %ylim([0.2 0.6]);
            pause
        end

    end

end