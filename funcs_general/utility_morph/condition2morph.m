


function [sch, morphA, morphB] = condition2morph(condition)

%     % Set up scheme A matrix (just take 2nd value in pair/unpaired)
%     SchAMat = [ 2 10 12 16 18  4 ;...
%                20 30 -1 -1 40 50 ;...
%                22 -1 32 42 -1 52 ;...
%                26 -1 46 36 -1 56 ;...
%                28 48 -1 -1 38 58 ;...
%                 6 60 62 66 68  8 ];
%     SchAMat(SchAMat == -1) = nan;
%     
%     % Set up Scheme B matrix
%     SchBMat = SchAMat + 68;
% 

    
    % Set up scheme A matrix (just take 2nd value in pair/unpaired)
    
    SchAMat(:,:,1)=[ 2 10 12 14 16 18  4 ;...
                    20 30 -1 -1 -1 40 50 ;...
                    22 -1 32 -1 42 -1 52 ;...
                    24 -1 -1 34 -1 -1 54 ;...
                    26 -1 46 -1 36 -1 56 ;...
                    28 48 -1 -1 -1 38 58 ;...
                     6 60 62 64 66 68  8 ];
                
    SchAMat(:,:,2)=[-1 -1 -1 -1 -1 -1 -1 ;...
                    -1 -1 -1 -1 -1 -1 -1 ;...
                    -1 -1 -1 -1 -1 -1 -1 ;...
                    -1 -1 -1 44 -1 -1 -1 ;...
                    -1 -1 -1 -1 -1 -1 -1 ;...
                    -1 -1 -1 -1 -1 -1 -1 ;...
                    -1 -1 -1 -1 -1 -1 -1 ];
    SchAMat(SchAMat == -1) = nan;
    
    % Set up Scheme B matrix
    SchBMat = SchAMat + 68;
                
%     
%     SchAMat
%     SchBMat
    
    % Round up to nearest odd number
    condition = ceil(condition / 2) * 2;
    
    % Calculate morph indices
    sz = size(SchAMat);
    if condition >= 2 && condition <= 68
        ind = find(SchAMat == condition);
        sch='A';
    elseif condition >= 70 && condition <= 136
        ind = find(SchBMat == condition);
        sch='B';
    else
        fprintf('Invalid condition number');
    end
    [I, J, ~] = ind2sub(sz,ind);
    
    [morphA] = convert_to_percent(I);
    [morphB] = convert_to_percent(J);
    
%     morphA = (7-I)*20;
%     morphB = (7-J)*20;

end


function morph = convert_to_percent (I)
    if I == 4; morph = 50;
    else
        if I >= 4; I=I-1; end
        morph = (6-I)*20;
    end
    
    
end



