   

function groups_merged = grouppairs_merge(group,iscirc,operation)
    % Does a diff on adjacent elements of group (i.e. subtracts them)
    %
    % Now replaced by grouppairs_diff
    
    if nargin < 3
        operation = [];      % 0 - Raw difference: ind1-ind2
                             % 1 - Percent difference: (ind1-ind2) / mean(ind2)
                             % 2 - Division for percent change: (ind1/mean(ind2) * 100)
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
    
    for i = 1:floor(length(group)/2)
        % Import default attributes
        groups_merged(i) = group(i*2-1);        % Inherit everything from 1st group in pair
        
        % % % % Do operation for data % % % % 
        if operation==0;
            z = group(i*2-1).data - group(i*2).data;
        elseif operation == 1
            warning('This takes % change in each trace individually; should not use this mode'); % This can result in singularities.
            sz = size(group(i*2).data);
            %z = (group(i*2-1).data - group(i*2).data)./ group(i*2).data * 100;
            z = (group(i*2-1).data - group(i*2).data)./ repmat(mean(group(i*2).data,2),[1,sz(2),1,1,1]) * 100;
            groups_merged(i).data_name = '% change';
        elseif operation == 2
            sz = size(group(i*2).data);
            z = group(i*2-1).data./ repmat(mean(group(i*2).data,2),[1,sz(2),1,1,1]) * 100;
            groups_merged(i).data_name = '% change';
            
        elseif operation == 4
                % Bootstrap normalized difference
                % E.g. Difference between data1 and data 2... 
                    % normalized by the MEAN of permuted distribution of the absolute difference.
                % Useful when group(1) and group(2) have different
                % numcells.
            data1 = group(2*i-1).data;
            data2 = group(2*i).data;
            
            my_function = @(x,y) abs(mean(x,2) - mean(y,2));
            data_permutednull = mean_of_permuted_null(my_function,data1,data2);

            z = (mean(data1,2)- mean(data2,2)) ./ data_permutednull;
            
            % % % % % % % % TEST! Force output to use null distribution  % % % % % % % % 
%             my_function2 = @(x,y) abs(mean(x,2) - mean(y,2));
%             Nboot = 1;
%             z = mean_of_permuted_null(my_function2,data1,data2,Nboot) ./ data_permutednull;
            % % % % % % % % End test  % % % % % % % % 
        end
            
        if iscirc
            ind=z>pi; z(ind)=z(ind)-2*pi;
            ind=z<-pi; z(ind)=z(ind)+2*pi;
            %z = abs(z);
        end
        groups_merged(i).data = z;
        
        
        % % % % Do operation for datastats % % % % 
        if operation==0;
            z = group(i*2-1).datastats - group(i*2).datastats;
        elseif operation == 1
            sz = size(group(i*2).datastats);
            %z = (group(i*2-1).datastats - group(i*2).datastats)./group(i*2).datastats * 100;
            z = (group(i*2-1).datastats - group(i*2).datastats)./repmat(mean(group(i*2).datastats,2),[1,sz(2),1,1,1]) * 100;
        elseif operation == 2
            sz = size(group(i*2).datastats);
            z = group(i*2-1).datastats./repmat(mean(group(i*2).datastats,2),[1,sz(2),1,1,1]) * 100;
            
        elseif operation == 4
                % Bootstrap normalized difference (via permutation test)
                % Useful when group(1) and group(2) have different
                % numcells.
                
            skip_datastats = 1;     % Set to 1 to skip calculating datastats (will be a single scalar anyways), in order to save time.
            if skip_datastats
                groups_merged(i).datastats = [];        
            else
            
                data1 = group(2*i-1).datastats;
                data2 = group(2*i).datastats;

                my_function = @(x,y) abs(mean(x,2) - mean(y,2));
                data_permutednull = mean_of_permuted_null(my_function,data1,data2);

                z = (mean(data1,2)- mean(data2,2)) ./ data_permutednull;
            end

            %group2(i).zlims_desired =[-2 2];
        end
        if iscirc
            ind=z>pi; z(ind)=z(ind)-2*pi;
            ind=z<-pi; z(ind)=z(ind)+2*pi;
            %z = abs(z);
        end
        groups_merged(i).datastats = z;
        
        
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