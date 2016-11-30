
function [t60, t60_irr, t50, t50_irr, t80, t80_irr, t100, t100_irr] = get_boundary_trials(md)

        Ntrials = length(md.start_times);

        only_along_active_boundary = 1;
        
        sch = md.sch;
        morphA = md.morphA;
        morphB = md.morphB;
        
        if only_along_active_boundary
            schA = strcmp(sch,'A');
            schB = strcmp(sch,'B');
        else
            schA = true(size(sch));
            schB = true(size(sch));
        end
        
        t50 = (schA & (morphA == 50)) | (schB & (morphB == 50));
        t50_irr = (schB & (morphA == 50)) | (schA & (morphB == 50));
        
        t60 = (schA & (morphA == 40 | morphA == 60)) | (schB & (morphB == 40 | morphB == 60));
        t60_irr = (schB & (morphA == 40 | morphA == 60)) | (schA & (morphB == 40 | morphB == 60));
        
        t80 = (schA & (morphA == 20 | morphA == 80)) | (schB & (morphB == 20 | morphB == 80));
        t80_irr = (schB & (morphA == 20 | morphA == 80)) | (schA & (morphB == 20 | morphB == 80));
        
        t100 = (schA & (morphA == 0 | morphA == 100)) | (schB & (morphB == 0 | morphB == 100));
        t100_irr = (schB & (morphA == 0 | morphA == 100)) | (schA & (morphB == 0 | morphB == 100));
        
        
%         if only_along_active_boundary
%             t60 = t60A | t60B;       % Uncomment for along active boundary
%             t60_irr = t60A_irr | t60B_irr;
%             t50 = t50A | t50B;
%             t50_irr = t50A_irr | t50B_irr;
%         else
%             t60 = (morphA >= 40 & morphA <=60) | (morphB >= 40 & morphB <=60);    % Uncomment for either boundary (SHOULDN'T THE MIDDLE AND BE AN OR?)
%         end
        
        % Reduce switch trials to be indexed by samples
        sot = md.sample_on_trials;
        t50 = t50(sot);
        t50_irr = t50_irr(sot);
        
        t60 = t60(sot);
        t60_irr = t60_irr(sot);
        
        t80 = t80(sot);
        t80_irr = t80_irr(sot);
        
        t100 = t100(sot);
        t100_irr = t100_irr(sot);
        
        
%         % Display percent switch trials
%         fprintf('Percent boundary trials = %g \n', (sum(tboundary) / length(tboundary) *100) );
end

