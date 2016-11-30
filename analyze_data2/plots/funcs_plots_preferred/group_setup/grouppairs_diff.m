   

function groups_merged = grouppairs_diff(group,iscirc,xlims_desired0)
    % Does a diff on adjacent elements of group (i.e. subtracts them)
    
    if nargin < 2
        iscirc = 0;
    end
    
    for i = 1:floor(length(group)/2)
        % Import default attributes
        groups_merged(i) = group(i*2-1);        % Inherit everything from 1st group in pair
        groups_merged(i).ctgs = floor(groups_merged(i).ctgs/2)+1;   % The CTGS will now match to appropariate entry in diffed pls.
        
        % Do subtraction of data
        if ~isempty(group(i*2-1).data) && ~isempty(group(i*2).data)
            z = group(i*2-1).data - group(i*2).data;
            if iscirc
                ind=z>pi; z(ind)=z(ind)-2*pi;
                ind=z<-pi; z(ind)=z(ind)+2*pi;
                %z = abs(z);
            end
            groups_merged(i).data = z;
        end
        
        
        % Do subtraction of datastats
        if ~isempty(group(i*2-1).datastats) && ~isempty(group(i*2).datastats)
            z = group(i*2-1).datastats - group(i*2).datastats;
            if iscirc
                ind=z>pi; z(ind)=z(ind)-2*pi;
                ind=z<-pi; z(ind)=z(ind)+2*pi;
                %z = abs(z);
            end
            groups_merged(i).datastats = z;
        end
        
        % Update legend (vs)
        n1 = group(i*2-1).legend;
        n2 = group(i*2).legend;
        
        % Trim n1 
        if ~isempty(n1)
            n1 = strrep(n1,'morphs','');
            n1 = strrep_if_in_n2(n1,'irr','',n2);
            n1 = strrep_if_in_n2(n1,'rel','',n2);
            n1 = strrep_if_in_n2(n1,'SchA','',n2);
            n1 = strrep_if_in_n2(n1,'SchB','',n2);
        end
        
        % Trim n2 RT
        if ~isempty(n2)
            n2 = strrep(n2,'Incongruent','');
            n2 = strrep(n2,'Congruent','');
            n2 = strrep(n2,'Ctg','');
        end
        
        n1n2 = [n1 ' vs ' n2];
        n1n2 = strrep(n1n2,'   ','  '); % Remove triple spaces
        n1n2 = strrep(n1n2,'  ',' '); % Remove double spaces
        groups_merged(i).legend = n1n2;
        
        
        % Update xlims_desired
        if exist('xlims_desired0','var'); groups_merged(i).xlims_desired = xlims_desired0; end
        
        
    end
    
    % Update data_name (ylabel)
        if ~isempty(group(1).data_name)
            groups_merged(1).data_name = ['\Delta ' groups_merged(1).data_name];
        end
    
end

function n1_out = strrep_if_in_n2(n1,str1,str2,n2)
    if ~isempty(strfind(n2,str1))
        n1_out = strrep(n1,str1,str2);
    else
        n1_out = n1;
    end
end