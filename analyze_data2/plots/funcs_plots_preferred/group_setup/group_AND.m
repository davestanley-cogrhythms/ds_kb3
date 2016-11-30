


function gr_out = group_AND(varargin)
    %gr_out = group_AND(gr1,gr2,gr3,...)
    % gr1...grN must be single structures, not arrays
    
    gr_arr = horzcat(varargin{:});

    gr_out = gr_arr(1);
    for i = 2:nargin
        gr_out = group_AND_pair (gr_out,gr_arr(i));
    end

end



function gr_out = group_AND_pair (gr1,gr2)
    % This "ANDs" gr2 criteria into gr
    % This essentially merges the non-2 criteria from both groups.
    % If the two groups have non-2 criteria at the same positions, displays
    % a warning.
    
    % All other fields of the output group are derived from Gr1.
    
    crit1=[gr1.criteria gr1.criteria_alt gr1.criteria_sfc];
    crit2=[gr2.criteria gr2.criteria_alt gr2.criteria_sfc];
    N1=size(crit1,1); N2=size(crit2,1); Ncrit = size(crit1,2);
    
    % Test to see if there is any discrepency between the criteria that are
    % the same
    overlap = 0;
    nonidentical_merge = 0;
    for i = 1:N1
        for j = 1:N2
            ccrit1 = crit1(i,:); ccrit2 = crit2(j,:);
            ind1=ccrit1 ~= 2; ind2 = ccrit2 ~=2;
            overlap_ind = ind1 & ind2;
            if any(overlap_ind);
                overlap = 1;
                fprintf('Overlap detected \n');
                if any(ccrit1(overlap_ind) ~= ccrit2(overlap_ind))
                    nonidentical_merge = 1;
                    fprintf(['Offending gr1 criteria: ' num2str(ccrit1) '\n']);
                    fprintf(['Offending gr2 criteria: ' num2str(ccrit2) '\n']);
                    warning('Overlapping portions of gr1 and gr2 criteria do not match');
                end
            end
        end
    end
    
    
    % Perform the "and"
    crit_merged = 2*ones(N1*N2,Ncrit);
    k=1;
    for i = 1:N1
        for j = 1:N2
            ccrit1 = crit1(i,:); ccrit2 = crit2(j,:);
            crit_merged(k,:) = ccrit2;
            ind1 = ccrit1 ~=2;
            crit_merged(k,ind1) = ccrit1(ind1);      % Move crit1 values into the merged criteria
            k=k+1;
        end
    end
    
    % Separate critera back into crit, crit_alt, and crit_sfc
    sz1=size(gr1.criteria); sz2=size(gr1.criteria_alt); sz3=size(gr1.criteria_sfc);
    critm_cell = mat2cell(crit_merged,[N1*N2],[sz1(2) sz2(2) sz3(2)]);
    critm = critm_cell{1}; critm_alt = critm_cell{2}; critm_sfc = critm_cell{3};
    
    gr_out(1:N2) = gr1;
    gr_out = group_merge(gr_out);
    for k = 1:N1*N2
        gr_out.criteria(k,:) = critm(k,:);
        gr_out.criteria_alt(k,:) = critm_alt(k,:);
        gr_out.criteria_sfc(k,:) = critm_sfc(k,:);
    end
    
    
    gr_out.legend = [gr1.legend ' & ' gr2.legend];
    
%     gr_arr = horzcat(varargin{:});
%     fieldn = fieldnames(gr_arr);
%     
%     for j = 1:length(fieldn)
%         gr_out.(fieldn{j}) = cat(1,gr_arr(:).(fieldn{j}));
%     end
end