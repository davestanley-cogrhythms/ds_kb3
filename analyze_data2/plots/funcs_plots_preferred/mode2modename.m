

function filename_out = mode2modename (mode)

    modes_name = num2str(mode,12);
    modes_name = [strrep(modes_name,'.','_')];
    filename_out = [modes_name];
end