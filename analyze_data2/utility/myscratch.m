


i=0;
i=i+1; sfc_mode_arr{i}=22.4018111;
% i=i+1; sfc_mode_arr{i}=22.4713111;
% i=i+1; sfc_mode_arr{i}=22.4718111;

parfor i = 1:length(sfc_mode_arr)
    %delete_unloadable_files(sfc_mode_arr{i},2);
    delete_unloadable_files(sfc_mode_arr{i},3);
end


if 0

    sfc_mode = 22.4415111; folder_prefix = 'test3'; delete_unloadable_files(sfc_mode,3,folder_prefix);
    sfc_mode = 22.4415111; delete_unloadable_files(sfc_mode,2,folder_prefix)
    sfc_mode = 22.4415111; folder_prefix = 'test4'; delete_unloadable_files(sfc_mode,3,folder_prefix);
    sfc_mode = 22.4415111; delete_unloadable_files(sfc_mode,2,folder_prefix)

end


%%

i=0;
i=i+1; sfc_mode_arr{i}=22.4018111;
i=i+1; sfc_mode_arr{i}=22.4713111;
i=i+1; sfc_mode_arr{i}=22.4718111;

parfor i = 1:length(sfc_mode_arr)
    count_i0j0_split_files(sfc_mode_arr{i},2);
    count_i0j0_split_files(sfc_mode_arr{i},3);
end


if 0

    sfc_mode = 22.4415111; folder_prefix = 'test3'; count_i0j0_split_files(sfc_mode,3,folder_prefix);
    sfc_mode = 22.4415111; folder_prefix = 'test3'; count_i0j0_split_files(sfc_mode,2,folder_prefix)
    sfc_mode = 22.4415111; folder_prefix = 'test4'; count_i0j0_split_files(sfc_mode,3,folder_prefix);
    sfc_mode = 22.4415111; folder_prefix = 'test4'; count_i0j0_split_files(sfc_mode,2,folder_prefix)

end