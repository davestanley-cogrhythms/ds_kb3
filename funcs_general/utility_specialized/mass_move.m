



x = dir;
src_dir_prefix = 'f';
src_filename = 'publish_run_analysis.pdf';
Npre=length(src_dir_prefix);
outprefix = 'SFC_';

for i = 1:length(x)
    if x(i).isdir && x(i).name(1:Npre) == src_dir_prefix
        [pathstr, name, ext] = fileparts(src_filename);
        mvsource = fullfile(x(i).name,src_filename);
        mvdest = [outprefix x(i).name(Npre+1:end) ext];
        movefile(mvsource,mvdest)
    end
end


