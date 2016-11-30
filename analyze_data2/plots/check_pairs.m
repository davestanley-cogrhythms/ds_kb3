

% Note: Need to run plots_preferred_unitpairs first!!
%% Check that the pairing of FFCs and Units is done correctly for random units

fn1 = wrkspc_buffer.sfc_2_20151.stage3.funames1D;
fn2 = wrkspc_buffer.sfc_22_40131.stage3.funames1D;
                        % NOTE: ALTHOUGH THE MYPAIRS WERE MISORDERED IN THE
                        % PREVIOUS DATA, THIS ACTUALLY DID NOT AFFECT MY
                        % RESULTS AT ALL. THE REASON FOR THIS IS BECAUSE
                        % EACH MYPAIR ELEMENT STILL CORRECTLY IDENTIFIED
                        % THE ELECTRODE NUMBER OF THE CORRESPONDING Cave
                        % varaible. The reason the funames1D were messed up
                        % is because funames1D is calculated separately, as
                        % and stored in the "firsts."
                        


random_pair = 1907;
random_pair = 579;
% random_pair = 18;
%random_pair = round(unifrnd(1,length(choup_ind)));  % Pick a random one of my 1976 FFC pairs (the ones associated with units)

names1 = [fn1{choup(random_pair,1)} ':' fn1{choup(random_pair,2)}]
names2 = fn2{choup_ind(random_pair)}
names3 = funames1D{random_pair}



%% Now for that random unit, check that md.unit_LFP_pairs is correct!

% First, get file number
fname = names1(1:7);
all_names = {md.recording_name};
ind = strfind(all_names,fname);
ind = ~cellfun(@isempty,ind);   % Index of current file

% Then, take the md corresponding to that file
mdc = md(ind);  % Current MD

% Calculate expected unit/LFP pairs
clear names
for i = 1:size(mdc.unit_LFP_pairs,2)
    names{i} = ['File:' fname 'Unit:' mdc.unit_names{mdc.unit_LFP_pairs(1,i)} 'vs LFP:' mdc.lfp_names{mdc.unit_LFP_pairs(2,i)}];
end

names'


%% Check that the pairing of FFC and Units is done correctly for all units

matches = [];
for i = 1:length(choup_ind)
    random_pair = i;
    
    % Get names of SFC and FFC electrode pairs
    sfc_names1 = fn1{choup(random_pair,1)};
    sfc_names2 = fn1{choup(random_pair,2)};
    sfc_names = [ sfc_names1 ':' sfc_names2];
    ffc_names = fn2{choup_ind(random_pair)};
    
    % Make sure file names are the same
    fname_match = strcmp(names1(1:11),names2(1:11));
    
    % Make sure electrode numbers are the same
        % First extract electrode numbers from SFC recordings  
    sfc_names;
    ind1 = strfind(sfc_names1,'PFC Elect') + length('PFC Elect');
    ind2 = strfind(sfc_names1,'unit') - 1;
    sfc_electnum1 = sfc_names1(ind1:ind2);
    sfc_electnum1 = str2num(sfc_electnum1);
    
    ind1 = strfind(sfc_names2,'PFC Elect') + length('PFC Elect');
    ind2 = strfind(sfc_names2,'unit') - 1;
    sfc_electnum2 = sfc_names2(ind1:ind2);
    sfc_electnum2 = str2num(sfc_electnum2);
    
        % Then extract electrode numbers from FFC recordings
    ffc_names;
    ind = strfind(ffc_names,':');
    ffc_names1 = ffc_names(1:ind-1);
    ffc_names2 = ffc_names(ind+1:end);
    
    ind = strfind(ffc_names1,'PFC Elect') + length('PFC Elect');
    ffc_electnum1 = ffc_names1(ind:end);
    ffc_electnum1 = str2num(ffc_electnum1);
    
    ind = strfind(ffc_names2,'PFC Elect') + length('PFC Elect');
    ffc_electnum2 = ffc_names2(ind:end);
    ffc_electnum2 = str2num(ffc_electnum2);
    
    electL_match = sfc_electnum1 == ffc_electnum1;
    electR_match = sfc_electnum2 == ffc_electnum2;
    
    matches = [matches; [fname_match electL_match electR_match]];
    
    
end

figure; plot(matches(:,2:3),'.');
figure; plot(choup_ind,matches(:,2:3),'.');

%% Check that the number of lfp_names in MD matches the number of channels in each file - everything looks good!

    path_lfp = fullfile(getpath('path_lfp'));
    [file_list] = get_filelist(path_lfp,'*.mat');
    
    nlfp_files = [];
    for i = 1:length(file_list)
        s=load(fullfile(path_lfp,file_list{i}));
        nlfp_files(i) = size(s.lfpmat_int16,2);
    end
    clear s
    
    nlfp_md = arrayfun(@(s) length(s.lfp_names),md);

    figure; plot([nlfp_files(:), nlfp_md(:)]);
    sum(nlfp_files(:) ~= nlfp_md(:))

    
   





