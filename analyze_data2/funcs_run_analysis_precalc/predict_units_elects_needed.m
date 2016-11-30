function [units_needed,elects_needed] = predict_units_elects_needed(ue_pairs,mypairs,j0)

    % Set up mypairs for iteration

    units_needed = [];
    elects_needed = [];
    switch ue_pairs
        case {0,2,3};
            
            
            units_needed=unique(mypairs(j0,1));
            elects_needed=unique(mypairs(j0,2));
            
            
        case 4                                                  % Current field versus all field (FFC)
            
            elects_needed=unique([mypairs(j0,1); mypairs(j0,2)]);
            
            
        case 5                                                  % Current unit versus all units (SSC)
            units_needed=unique([mypairs(j0,1); mypairs(j0,2)]);
        case 6                                                  % Current electrode only
            elects_needed=unique([mypairs(j0,1); mypairs(j0,2)]);
        case 7                                                  % Current unit only (non-parametric unit analysis)
            units_needed=unique([mypairs(j0,1); mypairs(j0,2)]);
    end


end