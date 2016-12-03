
classdef Grp

    properties
        
        % Minimal properties that must be specified (in instantiation)
        criteria
        criteria_alt
        criteria_sfc
        ctgs = [1];
        
        % Additional metadata
        coords
        right
        wrong
        diff
        score       % Arbitrary variable used to store some "score" or "ranking" value. Typically used for score2grouping
        netprefA
        metadata
        freqband_stats      % Frequency band used for calculating statistics
        timeband_stats      % Time band used for calculating statistics (if spectrogram)
        tf_label_stats
        freqband_perm       % Frequency band used for calculating permuted significance
        timeband_perm       % Time band used for calculating permuted significance (if spectrogram)
        tf_label_perm
        is_spectrogram      % Flag for spectrogram data
        
        % Legend settings
        legend                  % May contain latex
        legend_plaintext        % Plaintext alternative legend when specified interpreter is not available (some custom functions may need this)
        interpreter = 'none'   % Specifies interpreter to use
        legendtype = 'Cell Group';
        
        % Other metadata info
        data_type = 'FFC';              % Returned by mode2datatype function. 
        data_name = '';        % Ylabel
        xdata_name = 'Freq (Hz)';       % Xlabel
        xdata2_name = 'Time (s)';       % Xlabel2 (actually this is the y axis for spectrogram)
        
        Nwind = [];
        tapers = [];
        full_bandwidth = [];
        Nwind_perm = [];
        tapers_perm = [];
        full_bandwidth_perm = [];
        
        % Plotting stats
        xlims_desired = [0 120];
        ylims_desired
        zlims_desired
        colour_desired = [];
        colourmap_desired = [];
        
        % Populated data
        cells
        data
        data_STE
        data_pvals
        xdata
        xdata2
        datastats
        datastats_STE
        datastats_pvals
        data_overlay1       % For spectrograms, this will be transparency matrix
        data_overlay2       % For spectrograms, this will be contours matrix
        numcells = -1
        funames
    end
    
    methods
        
        function gr = Grp
            N_ctgs_base = get_Nctgs_base;
            N_ctgs_extra = get_Nctgs_extra;
            N_ctgs = N_ctgs_base + N_ctgs_extra;
            N_pairs = N_ctgs / 2;

            gr.criteria = 2*ones(1,N_pairs);
            gr.criteria_alt = gr.criteria;
            gr.criteria_sfc = [2, 2];   % First entry is for SFC, second is for coding for animal (cheap hack =/)
            gr.ctgs = [1];
        end
        
        function gr_sd = Grp_SD(gr)     % Converter method; instructions to convert to Subclass Grp_SD
            gr_sd = Grp_SD;            
            gr_sd = inheritObj(gr_sd,gr);
        end
        
        function gr = reset(gr)
            fprintf('test \n');
        end
        
        
        function gr = append_string2(gr,myfield,mytext)
            %Adds text to an existing string in object

            gr.(myfield) = [gr.(myfield) mytext];

        end
        
        function gr = prepend_string(gr,myfield,mytext)
            %Adds text to an existing string in object

            gr.(myfield) = [mytext,gr.(myfield) ];

        end

        
    end

end


function output = inheritObj(output,input)
    % Merges contents of input into output.
    C = metaclass(input);
    P = C.Properties;
    for k = 1:length(P)
        if ~P{k}.Dependent
            output.(P{k}.Name) = input.(P{k}.Name);
        end
    end
end


