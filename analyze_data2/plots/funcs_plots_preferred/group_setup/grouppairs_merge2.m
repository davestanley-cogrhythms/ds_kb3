   

function groups_merged = grouppairs_merge2(group,iscirc,operation)
%% function groups_merged = grouppairs_merge2(group,iscirc,operation)
    % Does a diff on adjacent elements of group (i.e. subtracts them)
    %
    % Now replaced by grouppairs_diff
    
    if nargin < 3
        operation = [];      % 0 - Raw difference: ind1-ind2
                             % 1 - Percent difference: (ind1-ind2) / mean(ind2)
                             % 2 - Percent difference (ind1-ind2) / mean((ind1+ind2)/2)
                             % 4 - Difference divided by z-score
    end
    
    if nargin < 2
        iscirc = [];
    end
    
    if isempty(operation)
        operation = 0;
    end
    
    if isempty(iscirc)
        iscirc = 0;
    end
    
    myalpha = 0.01;
    
    for i = 1:floor(length(group)/2)
        % Import default attributes
        groups_merged(i) = group(i*2-1);        % Inherit everything from 1st group in pair
        
        % Perform operation
        group1 = group(2*i-1); group2 = group(2*i);
        ispaired = get_ispaired(group1,group2);
        if ispaired
            [groups_merged(i).data] = do_operation_paired(group1,group2,'data',operation,iscirc);
            [groups_merged(i).datastats] = do_operation_paired(group1,group2,'datastats',operation,iscirc);  % should probably just get rid of this.
        else
            [z_mu, z_ste, pvals] = do_operation_unpaired(group1,group2,'data',operation,iscirc);
            groups_merged(i).data_mu = z_mu;
            groups_merged(i).data_STE = z_ste;
            groups_merged(i).data_pvals = pvals;
            [groups_merged(i).datastats] = do_operation_unpaired(group1,group2,'datastats',operation,iscirc);  % should probably just get rid of this.
        end

        % Update legend (vs)
        n1 = group(i*2-1).legend;
        n2 = group(i*2).legend;
        
        if ~isvector(n1); n1 = n1(1,:); end
        if ~isvector(n2); n2 = n2(1,:); end
        
        % Trim n1 
        n1 = strrep(n1,'morphs','');
        n1 = strrep(n1,'irr.','');
        n1 = strrep(n1,'rel.','');
        
        % Trim n2 RT
        n2 = strrep(n2,'Incongruent','');
        n2 = strrep(n2,'Congruent','');
        
        % Merge
        nout = [n1 ' vs ' n2];
        
        % Remove double/triple spaces
        nout = strrep(nout,'   ',' ');
        nout = strrep(nout,'  ',' '); 
        
        groups_merged(i).legend = [nout];
        
        
    end
    
end


function data_permutednull = mean_of_permuted_null(my_function,data1,data2,Nboot)

    if nargin < 4
        Nboot = 100;
    end

    Ncells1 = size(data1,2);
    Ncells2 = size(data2,2);


    data_pooled = cat(2,data1,data2);
    data_permutednull = zeros(size(data_pooled(:,1,:,:,:,:,:)));
    for i = 1:Nboot
        currperm = randperm(Ncells1+Ncells2);
        data_perm1 = data_pooled(:,currperm(1:Ncells1),:,:);        % Draw random samples from pooled data for data1
        data_perm2 = data_pooled(:,currperm(Ncells1+1:end),:,:);    % Draw random samples from pooled data for data2

        % Apply our function to the data and then add it to our
        % running tally. In order to conserve space, we keep a
        % running tally rather than saving every single
        % permutation.
        data_permutednull = data_permutednull + my_function(data_perm1,data_perm2);   % this is summing things up; we will later take the mean
    end


    data_permutednull = data_permutednull / Nboot;      % Finally, divide the sum so as to get the average
end


function [z, z_ste, pvals] = do_operation_unpaired(group1,group2, fieldname, operation,iscirc)
%% function [z, z_ste, z_pval] = do_operation_unpaired(group1,group2, fieldname, operation, iscirc)

        z = [];
        z_ste = [];
        pvals = [];
        
        data1 = group1.(fieldname);
        data2 = group2.(fieldname);
        
        switch operation
            case 0
                % Mean difference
                z = mean(data1,2) - mean(data2,2);

                % Calculate standard error (assume gaussian!)
                z_ste = calc_ste(data1,data2);

                % Calculate p values
                statsfunc = @ranksum;   % Since it's unpaired!
                pvals = get_pvals(data1,data2,statsfunc);

            case 1
                warning('Untested!')
                % Mean % difference
                z = (mean(data1,2) - mean(data2,2)) ./ mean(data2,2) * 100;

                % Calculate standard error (assume gaussian!)
                z_ste = calc_ste(data1,data2) ./ mean(data2,2) * 100;

                % Calculate p values
                statsfunc = @ranksum;   % Since it's unpaired!
                pvals = get_pvals(data1,data2,statsfunc);

            case 2
                warning('Untested!')
                % Mean % difference
                z = (mean(data1,2) - mean(data2,2)) ./ (mean(data2,2) + mean(data2,2))/2 * 100;

                % Calculate standard error (assume gaussian!)
                z_ste = calc_ste(data1,data2) ./ (mean(data2,2) + mean(data2,2))/2 * 100;

                % Calculate p values
                statsfunc = @ranksum;   % Since it's unpaired!
                pvals = get_pvals(data1,data2,statsfunc);

            case 4
                % Bootstrap normalized difference
                % E.g. Difference between data1 and data 2... 
                    % normalized by the MEAN of permuted distribution of the absolute difference.
                % Useful when group(1) and group(2) have different
                % numcells.

                my_function = @(x,y) abs(mean(x,2) - mean(y,2));
                data_permutednull = mean_of_permuted_null(my_function,data1,data2);

                z = (mean(data1,2)- mean(data2,2)) ./ data_permutednull;

                % Calculate standard error (assume gaussian!)
                z_ste = calc_ste(data1,data2) ./ data_permutednull;

                % Calculate p values
                statsfunc = @ranksum;   % Since it's unpaired!
                pvals = get_pvals(data1,data2,statsfunc);

            otherwise
                error('Undefined operation');

        end
        
        if iscirc && any(operation == [1,2,3])
            ind=z>pi; z(ind)=z(ind)-2*pi;
            ind=z<-pi; z(ind)=z(ind)+2*pi;
            %z = abs(z);
        end
end

function [z] = do_operation_paired(group1,group2, fieldname, operation,iscirc)
%% function [z, z_ste, z_pval] = do_operation_paired(group1,group2, fieldname, operation,iscirc)

        z = [];
        
        data1 = group1.(fieldname);
        data2 = group2.(fieldname);
        
         switch operation
            case 0
                z = data1 - data2;
            case 1
                sz = size(data2);
                %z = (data1 - data2)./ data2 * 100;
                z = (data1 - data2)./ repmat(mean(data2,2),[1,sz(2),1,1,1]) * 100;
            case 2
                sz = size(data2);
                z = (data1 - data2)./ repmat(mean((data1 + data2)/2,2),[1,sz(2),1,1,1]) * 100;
            otherwise
                error('Undefined operation');
        end

        if iscirc && any(operation == [1,2,3])
            ind=z>pi; z(ind)=z(ind)-2*pi;
            ind=z<-pi; z(ind)=z(ind)+2*pi;
            %z = abs(z);
        end

end


function pvals = get_pvals(data1,data2,statsfunc)
%% function pvals = get_pvals(data1,data2,statsfunc)
    if ndims(data1) == 2    % mat2cell can't handle dimension mismatches well.
        datac1 = mat2cell(data1,ones(1,size(data1,1)),size(data1,2));
        datac2 = mat2cell(data2,ones(1,size(data2,1)),size(data2,2));
    else
        datac1 = mat2cell(data1,ones(1,size(data1,1)),size(data1,2),ones(1,size(data1,3)));
        datac2 = mat2cell(data2,ones(1,size(data2,1)),size(data2,2),ones(1,size(data2,3)));
    end
    pvals = cellfun(statsfunc,datac1,datac2);
end


function z_ste = calc_ste(data1,data2)
    varmean1 = var(data1,[],2) / size(data1,2);  % Variance of mean of data1
    varmean2 = var(data2,[],2) / size(data2,2);  % Variance of mean of data2
    z_ste = sqrt(varmean1 + varmean2);           % Variance of data1-data2
end

function ispaired = get_ispaired(group1,group2)
    % Check if data is paired or unpaired
    if all(group1.criteria(:) == group2.criteria(:)) 
        ispaired = true;
    else
        ispaired = false;
    end
end