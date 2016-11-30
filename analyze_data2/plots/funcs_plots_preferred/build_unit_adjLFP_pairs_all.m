

function unit_adjLFP = build_unit_adjLFP_pairs_all(md)
    % md - array of md structures (generally spanned over file range, i.e. 1:79)

    %%
    Nunits = arrayfun(@(x) length(x.unit_names),md);
    Nelects = arrayfun(@(x) length(x.lfp_names),md);
    
    units = arrayfunu(@(n) 1:n,Nunits);
    
    func2 = @(m,x) arrayfun(@(x) get_adjacent_electrode(m,x),x);
    unit_adjLFP = cellfunu(func2, num2cell(md),units);
    
    % Equiv for loop
%     clear unit_adjLFP2
%     for i = 1:length(md)
%         currarr = [];
%         for j = 1:length(units{i})
%             currarr(j)=get_adjacent_electrode(md(i),units{i}(j));
%         end
%         unit_adjLFP2{i} = currarr;
%     end
    %%
    %unit_adjLFP2 = horzcat(unit_adjLFP2{:});
    
    
    unit_adjLFP = horzcat(unit_adjLFP{:});
    
    %%
    
    %Nunitsc=cumsum([0 Nunits(1:end-1)]);
    Nelectsc=cumsum([0 Nelects(1:end-1)]);
    
    %Nunitsc = cellfunu(@(x,y) repmat(x,[1,y]),num2cell(Nunitsc),num2cell(Nunits));
    Nelectsc = cellfunu(@(x,y) repmat(x,[1,y]),num2cell(Nelectsc),num2cell(Nunits));
    
    %Nunitsc = horzcat(Nunitsc{:});
    Nelectsc = horzcat(Nelectsc{:});

    
    unit_adjLFP = unit_adjLFP + Nelectsc;
    
%     
%     unit_LFP_pairs_all = horzcat(md.unit_LFP_pairs);
%     unit_LFP_pairs_all(1,:) = unit_LFP_pairs_all(1,:) + Nunitsc;
%     unit_LFP_pairs_all(2,:) = unit_LFP_pairs_all(2,:) + Nelectsc;
    
    % still to do - debug, comment, and remove for loop

end