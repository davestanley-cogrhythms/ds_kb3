

function outname = stagename(curr_stage)


    outname = num2str(curr_stage);
    outname = strrep(outname,'-','m');
    outname = ['stage' outname];

end