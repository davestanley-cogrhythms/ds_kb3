

function group_scored = score2grouping(group,score,order)
    
    if nargin < 3
        order = 'descend';
    end

    us = sort(unique(score),order); %unique scores
%     Nscoredgroups = length(us);
%     Ncriteria = size(group(1).criteria,1);
%     Nctgs = size(group(1).ctgs,1);
    
%     group_scored=struct('criteria',[zeros(Nscoredgroups,Ncriteria)],'ctgs',[zeros(Nscoredgroups,Nctgs)],'right',0,'wrong',0,'coords',{0});
%     % This doesn't work because the size of these zero arrays varies
%     group_scored = repmat(1,Nscoredgroups);
    

    for i = 1:length(us)
        index = find(us(i) == score);
        
        group_scored(i) = group_OR(group(index));
        group_scored(i).legend=['Pr ' num2str(us(i))];
    end
end


