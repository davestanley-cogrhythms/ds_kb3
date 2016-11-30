

function table_out = getcrit(group)
    % Dumps out a cell array
    % Each cell array entry contains the value of group.field at the
    % corresponding value of group.coords
    
    if nargin < 3; print_on = 1; end
    
    field = 'criteria';
    table_out1 = vertcat(group(:).(field));
    field = 'criteria_alt';
    table_out2 = vertcat(group(:).(field));
    field = 'criteria_sfc';
    table_out3 = vertcat(group(:).(field));
    
    table_out = [table_out1 table_out2 table_out3];
    
    
%     if print_on
%         table_out
%     end
    
end
