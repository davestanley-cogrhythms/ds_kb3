
% Will use all the for loops I want, since it's a short code
% This function loads Jefferson's preferred unit categories and stores them
% in the metadata array

function [md_struct, md_matrix] = calc_metadata_preferred(filename,md,p)

if ~exist('p','var')
    preferred_path = fullfile(getpath('path_preferred'),'PFC_category_unit_names_CRC.mat');
    p = load(preferred_path);    
end


unit_groupings = fieldnames(p);

Nunits = length(md.unit_names);
Ngroups = length(unit_groupings);

md_struct = build_md_expanded(filename,md);
md_matrix = zeros(Nunits,Ngroups);
for j = 1:Ngroups
    out_group = convert_grouping(getfield(p,unit_groupings{j}));
    [md_struct md_matrix(:,j)] = addgroup(md_struct,out_group,unit_groupings{j}); 
end



end

function md2 = build_md_expanded(filename,md)

    [pathstr, name, ext] = fileparts(filename);
    Nunits = length(md.unit_names);

    for i = 1:Nunits
        md2{i}.filename = name;
        
        electname = md.metadata.unit_LFP_pairs_metadata{i}.LFP_str;
        electind = strfind(electname,'_Elect_') + length('_Elect_');
        md2{i}.electnum = str2num(electname(electind:end));
        
        unitname = md.metadata.unit_LFP_pairs_metadata{i}.unit_str;
        unitind = strfind(unitname,'_unit_') + length('_unit_');
        md2{i}.unitnum = str2num(unitname(unitind:end));
    end

end


function out_group = convert_grouping(in_group)
    
    Nunits2 = size(in_group,1);
    
    for i = 1:Nunits2
        indelect = strfind(in_group(i,:),'_Elect_');
        indunit = strfind(in_group(i,:),'_unit_');
        
        ind2 = indelect-1;
        out_group{i}.filename = in_group(i,1:ind2);
        
        ind1 = indelect + length('_Elect_');
        ind2 = indunit - 1;
        out_group{i}.electnum = str2num(in_group(i,ind1:ind2));
        
        ind1 = indunit + length('_unit_');
        out_group{i}.unitnum = str2num(in_group(i,ind1:end));
    end
end

function [md2 mdarray] = addgroup(md2,grouping,groupname)
   Nunits = length(md2);
   Nunits_in_group = length(grouping);
   
   mdarray = zeros(Nunits,1);
   
   for i = 1:Nunits
       md2{i} = setfield(md2{i},groupname,0);
       for j = 1:Nunits_in_group
%            i
%            j
%            md2{i}.electnum
%            grouping{j}.electnum 
           if strcmp(md2{i}.filename,grouping{j}.filename) && md2{i}.electnum == grouping{j}.electnum && md2{i}.unitnum == grouping{j}.unitnum
               md2{i} = setfield(md2{i},groupname,1);
               mdarray(i,1) = 1;
           end
       end
   end
end




