

function unit_LFP_pairs_all = build_unit_LFP_pairs_all(md)
    % md - array of md structures (generally spanned over file range, i.e. 1:79)

    %%
    
    Nunits = arrayfun(@(x) length(x.unit_names),md);
    Nelects = arrayfun(@(x) length(x.lfp_names),md);
    Nunitsc=cumsum([0 Nunits(1:end-1)]);
    Nelectsc=cumsum([0 Nelects(1:end-1)]);
    
    Nunitsc = cellfunu(@(x,y) repmat(x,[1,y]),num2cell(Nunitsc),num2cell(Nunits));
    Nelectsc = cellfunu(@(x,y) repmat(x,[1,y]),num2cell(Nelectsc),num2cell(Nunits));
    
    Nunitsc = horzcat(Nunitsc{:});
    Nelectsc = horzcat(Nelectsc{:});
    
    unit_LFP_pairs_all = horzcat(md.unit_LFP_pairs);
    unit_LFP_pairs_all(1,:) = unit_LFP_pairs_all(1,:) + Nunitsc;
    unit_LFP_pairs_all(2,:) = unit_LFP_pairs_all(2,:) + Nelectsc;
    

end