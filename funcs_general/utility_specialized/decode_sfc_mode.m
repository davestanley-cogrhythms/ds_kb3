
function [mode_group, mode_subgroup] = decode_sfc_mode(sfc_mode)  % Returns pre-decimal and post-decimal values in sfc_mode (digits as array elements)
    myprecision = 12; % If go over size of double, we'll get random numbers
    y = num2str(sfc_mode,myprecision);
    z = arrayfunu(@(x) str2num(x),y);
    dec_ind = strfind(y,'.');
    mode_group = cell2mat(z(1:dec_ind-1));
    mode_subgroup_temp = cell2mat(z(dec_ind:end));
    
    mode_subgroup = zeros(1,myprecision);
    mode_subgroup(1,1:length(mode_subgroup_temp)) = mode_subgroup_temp;
    
end

