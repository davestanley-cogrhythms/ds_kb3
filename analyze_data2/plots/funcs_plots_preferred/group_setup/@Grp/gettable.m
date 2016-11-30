

function table_out = gettable(group, field,print_on)
    % Dumps out a cell array
    % Each cell array entry contains the value of group.field at the
    % corresponding value of group.coords
    
    if nargin < 3; print_on = 1; end
    
    allcoords = {group(:).coords};
    allcoords = vertcat(allcoords{:});
    acu = unique(allcoords,'rows');
    imax = max(acu(:,1));
    jmax = max(acu(:,2));
    
    if size(allcoords,1) ~= size(acu,1); contains_duplicates = true; 
    else contains_duplicates = false;end

    %if ~contains_duplicates
    if length(group) > 3 && ~contains_duplicates        % If there are no duplicates, we can merge into 1 table. Only do this if there would otherwise be more than 2 possible tables
        table_out = cell(imax,jmax);
        for i = 1:length(group)
            
            table_out = calctable_single(group(i),field,table_out);
        end
        table_out = {table_out};
    else
        for i = 1:length(group)
            table_out{i} = cell(imax,jmax);
            table_out{i} = calctable_single(group(i),field,table_out{i});
        end
    end
    
    
    if print_on
        fprintf('------------ \n');
        for l = 1:length(table_out)
            titlestr = ['Group ' num2str(l) ' ' group(l).legend '\n'];
            titlestr = strrep(titlestr,'\cup',' ');
            fprintf(titlestr);
            fprintf('------------ \n');
            currtab = table_out{l};
            for i = 1:size(currtab,1)
                for j = 1:size(currtab,2)
                    for k = 1:length(currtab{i,j})
                        fprintf(num2str(currtab{i,j}(k)));
                    end
                    fprintf('\t');
                end
                fprintf('\n');
            end
            fprintf('------------ \n');
        end
    end
end


function tab = calctable_single(gr,field,tab)
    
    grc = gr.coords;
    w = gr.(field);
    if size(grc,1) ~= size(w,1)
        w = w';                 % Make sure we're running along the correct dimension
    end
    for j = 1:size(grc,1)
        
        tab{grc(j,1),grc(j,2)} = w(j,:);
    end
end


