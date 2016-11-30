

function unit_LFP_pairs = get_allunit_LFP_mapping(md)
    
    % Build 1Du matrix (i.e. one entry for each of 549 units) of file numbers and unit numbers
    units_per_file = arrayfun(@(x) length(x.unit_names),md);
    tot_units =  sum(units_per_file);
    fnums_1Du = arrayfunu(@(x,y) repmat(x,1,y),1:length(md),units_per_file);
    fnums_1Du = cell2mat(fnums_1Du);
    unums_1Du = arrayfunu(@(x) 1:x, units_per_file); unums_1Du = cell2mat(unums_1Du);
    
    % Build 1Du matrix of unit numbers and lfp numbers
    unums_1Du_from_pairs = arrayfunu(@(x) x.unit_LFP_pairs(1,:),md); unums_1Du_from_pairs = cell2mat(unums_1Du_from_pairs); % unit numbers
    lnums_1Du_from_pairs = arrayfunu(@(x) x.unit_LFP_pairs(2,:),md); lnums_1Du_from_pairs = cell2mat(lnums_1Du_from_pairs); % lfp numbers
    if sum(unums_1Du ~= unums_1Du_from_pairs) > 0 fprintf('Sanity check - these two vectors should always be the same. Since they are not, something is wrong. \n'); end
    
    % % Do the same thing for lfp
    % % I.e. build 1Dl matrix (i.e. one entry for each of 873 lfp's) of file numbers and lfp numbers
    lfp_per_file = arrayfun(@(x) length(x.lfp_names),md);
    tot_lfp =  sum(lfp_per_file);
    fnums_1Dl = arrayfunu(@(x,y) repmat(x,1,y),1:length(md),lfp_per_file);
    fnums_1Dl = cell2mat(fnums_1Dl);
    lnums_1Dl = arrayfunu(@(x) 1:x, lfp_per_file); lnums_1Dl = cell2mat(lnums_1Dl);
    
    
    % % Pair them up
    my1Du = [fnums_1Du; unums_1Du; lnums_1Du_from_pairs];
    my1Dl = [fnums_1Dl; lnums_1Dl];

    % % Get the 1Dl index for each unit
    index = return_1Dl_index(my1Du(:,1), my1Dl); %Test
    u_indices = 1:tot_units;
    l_indices = arrayfun(@(x) return_1Dl_index(my1Du(:,x),my1Dl),u_indices);
    
    unit_LFP_pairs = [u_indices; l_indices];
    
end


function index = return_1Dl_index(curr_1Du, my1Df)

    index1 = curr_1Du(1) == my1Df(1,:);
    index2 = curr_1Du(3) == my1Df(2,:);
    index = [index1; index2];
    index = all(index,1);
    index = find(index);

end

