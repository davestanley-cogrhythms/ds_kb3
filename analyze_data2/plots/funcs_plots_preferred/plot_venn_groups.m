

function plot_venn_groups(group)

%          For a three circle diagram:
%            data is a seven element vector of:
%               |A|
%               |A and B|
%               |B|
%               |B and C|
%               |C|
%               |C and A|
%               |A and B and C|
%

work_in_percents = 0;
fixed_sizes = 1;

if length(group) > 2
    A = any(group(1).cells,2);
    B = any(group(2).cells,2);
    C = any(group(3).cells,2);
    
    mylegend = {group(1:3).legend};
    plot_venn(A,B,C,mylegend,work_in_percents,fixed_sizes);
else
    A = any(group(1).cells,2);
    B = any(group(2).cells,2);
    
    mylegend = {group(1:2).legend};
    plot_venn(A,B,mylegend,work_in_percents,fixed_sizes);
end


end


