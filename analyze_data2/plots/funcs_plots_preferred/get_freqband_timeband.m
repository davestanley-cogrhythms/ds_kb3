

function [s] = get_freqband_timeband(sfc_mode,opts_exclude)


    if nargin < 2
        animal_mode=0;  % Animal_mode: 0-both; 1-Animal L; 2-Animal O
    end
    
    i=1;
    s(i).timeband = [.6];
    s(i).freqband = [18];
    s(i).label = 'Default';
    s(i).label_short = 'AAA';
    
    i=0;
    if any(sfc_mode == [45.6018111011])
        if ~opts_exclude.excludeL
            % % % % % % % % % % % % % % % % % % % % Scheme A % % % % % % % % % % % % % % % % % % % % 
            i=i+1;
                s(i).label_short = 'LSA1';
                s(i).timeband = [.451];  % 
                s(i).freqband = [68.36];
                s(i).label = 'Animal L Sample Cat Gamma';
                
            i=i+1;
                s(i).label_short = 'LSA1';
                s(i).timeband = [.421];  % 
                s(i).freqband = [25.39];
                s(i).label = 'Animal L Sample Dog HighBeta';

            i=i+1;
                s(i).label_short = 'LDA1';
                s(i).timeband = [1.081];  
                s(i).freqband = [12.7];
                s(i).label = 'Animal L Delay Cat Alpha';

            i=i+1;
                s(i).label_short = 'LDA2';
                s(i).timeband = [.901];  % 
                s(i).freqband = [27.34];
                s(i).label = 'Animal L Delay Cat (small red dot) HighBeta';

            % % % % % % % % % % % % % % % % % % % % Scheme B % % % % % % % % % % % % % % % % % % % % 
            i=i+1;
                s(i).label_short = 'LSB1';
                s(i).timeband = [.481];  % 
                s(i).freqband = [24.41];
                s(i).label = 'Animal L Sample Tad Beta';

            i=i+1;
                s(i).label_short = 'LSB2';
                s(i).timeband = [.481];  % 
                s(i).freqband = [14.65];
                s(i).label = 'Animal L Sample Tad Alpha';
        end
        if ~opts_exclude.excludeO
            % % % % % % % % % % % % % % % % % % % % Scheme A % % % % % % % % % % % % % % % % % % % % 
            i=i+1;
                s(i).label_short = 'OSA1';
                s(i).timeband = [.301];   % 
                s(i).freqband = [15.62];
                s(i).label = 'Animal O Sample Dog Beta';

            i=i+1; 
                s(i).label_short = 'OSA2';
                s(i).timeband = [.181];   % 
                s(i).freqband = [7.812];
                s(i).label = 'Animal O Sample Cat Alpha';

            % % % % % % % % % % % % % % % % % % % % Scheme B % % % % % % % % % % % % % % % % % % % % 
            i=i+1;
                s(i).label_short = 'OSB1';
                s(i).timeband = [.391];   % 
                s(i).freqband = [25.39];
                s(i).label = 'Animal O Sample Goc HighBeta';
            
            i=i+1;
                s(i).label_short = 'OSB2';
                s(i).timeband = [.481];   % 
                s(i).freqband = [9.766];
                s(i).label = 'Animal O Sample Tad Alpha';

            i=i+1;
                s(i).label_short = 'ODB1';
                s(i).timeband = [.781];   % 
                s(i).freqband = [5.859];
                s(i).label = 'Animal O Delay Tad Theta';

            i=i+1;
                s(i).label_short = 'ODB2';
                s(i).timeband = [1.441];   % 
                s(i).freqband = [12.7];
                s(i).label = 'Animal O Delay Goc Alpha';

            i=i+1;
                s(i).label_short = 'ODB3';
                s(i).timeband = [1.381];   % 
                s(i).freqband = [23.44];
                s(i).label = 'Animal O Delay Goc HighBeta';

            i=i+1;
                s(i).label_short = 'ODB4';
                s(i).timeband = [1.201];   % 
                s(i).freqband = [43.95];
                s(i).label = 'Animal O Delay Goc Gamma';
        end
        
        
    elseif any(sfc_mode == [45.6013111010, 45.6013111011,45.6013111013])
        if ~opts_exclude.excludeL
            % % % % % % % % % % % % % % % % % % % % Scheme A % % % % % % % % % % % % % % % % % % % % 
            i=i+1;
                s(i).label_short = 'LSA1';
                s(i).timeband = [.451];  % 
                s(i).freqband = [68.36];
                s(i).label = 'Animal L Sample Cat Gamma';
                
            i=i+1;
                s(i).label_short = 'LSA1';
                s(i).timeband = [.481];  % 
                s(i).freqband = [25.39];
                s(i).label = 'Animal L Sample Dog HighBeta';

            i=i+1;
                s(i).label_short = 'LDA1';
                s(i).timeband = [.931];  
                s(i).freqband = [13.67];
                s(i).label = 'Animal L Delay Cat Alpha';

            i=i+1;
                s(i).label_short = 'LDA2';
                s(i).timeband = [.781];  % 
                s(i).freqband = [27.34];
                s(i).label = 'Animal L Delay Cat (small red dot) HighBeta';

            % % % % % % % % % % % % % % % % % % % % Scheme B % % % % % % % % % % % % % % % % % % % % 
            i=i+1;
                s(i).label_short = 'LSB1';
                s(i).timeband = [.511];  % 
                s(i).freqband = [25.39];
                s(i).label = 'Animal L Sample Tad Beta';

            i=i+1;
                s(i).label_short = 'LSB2';
                s(i).timeband = [.511];  % 
                s(i).freqband = [15.62];
                s(i).label = 'Animal L Sample Tad Alpha';
        end
        if ~opts_exclude.excludeO
            % % % % % % % % % % % % % % % % % % % % Scheme A % % % % % % % % % % % % % % % % % % % % 
            i=i+1;
                s(i).label_short = 'OSA1';
                s(i).timeband = [.331];   % 
                s(i).freqband = [17.58];
                s(i).label = 'Animal O Sample Dog Beta';

            i=i+1; 
                s(i).label_short = 'OSA2';
                s(i).timeband = [.301];   % 
                s(i).freqband = [5.859];
                s(i).label = 'Animal O Sample Cat Theta';

            % % % % % % % % % % % % % % % % % % % % Scheme B % % % % % % % % % % % % % % % % % % % % 
            
            i=i+1;
                s(i).label_short = 'OSB1';
                s(i).timeband = [.391];   % 
                s(i).freqband = [25.39];
                s(i).label = 'Animal O Sample Goc HighBeta';
                
            i=i+1;
                s(i).label_short = 'OSB2';
                s(i).timeband = [.481];   % 
                s(i).freqband = [9.766];
                s(i).label = 'Animal O Sample Tad Alpha';

            i=i+1;
                s(i).label_short = 'ODB1';
                s(i).timeband = [.721];   % 
                s(i).freqband = [5.859];
                s(i).label = 'Animal O Delay Tad Theta';

            i=i+1;
                s(i).label_short = 'ODB2';
                s(i).timeband = [1.441];   % 
                s(i).freqband = [11.72];
                s(i).label = 'Animal O Delay Goc Alpha';

            i=i+1;
                s(i).label_short = 'ODB3';
                s(i).timeband = [1.321];   % 
                s(i).freqband = [23.44];
                s(i).label = 'Animal O Delay Goc HighBeta';

            i=i+1;
                s(i).label_short = 'ODB4';
                s(i).timeband = [.991];   % 
                s(i).freqband = [46.88];
                s(i).label = 'Animal O Delay Goc Gamma';
                
            i=i+1;
                s(i).label_short = 'ODB5';
                s(i).timeband = [1.171];   % 
                s(i).freqband = [17.58];
                s(i).label = 'Animal O Delay Tad Beta';
        end
        
        
% % % % % % % % % % % % % % % % % % % % % % % % %         ////////
        

    elseif any(sfc_mode == [23.4013111010, 23.4013111011, 23.4013111013, 23.4018111013, 23.4018111011])
        if ~opts_exclude.excludeL
            % % % % % % % % % % % % % % % % % % % % Scheme A % % % % % % % % % % % % % % % % % % % % 
            i=i+1;
                s(i).label_short = 'LSA1';
                s(i).timeband = [.391];  % 
                s(i).freqband = [19.53];
                s(i).label = 'Animal L Sample Cat Beta';
                
            i=i+1;
                s(i).label_short = 'LSA2';
                s(i).timeband = [.481];  % 
                s(i).freqband = [1.953];
                s(i).label = 'Animal L Sample Cat Delta';

            i=i+1;
                s(i).label_short = 'LDA1';
                s(i).timeband = [1.171];  
                s(i).freqband = [19.53];
                s(i).label = 'Animal L Delay Dog Beta';

            i=i+1;
                s(i).label_short = 'LDA2';
                s(i).timeband = [.841];  % 
                s(i).freqband = [13.67];
                s(i).label = 'Animal L Delay Cat Alpha';
                
            i=i+1;
                s(i).label_short = 'LDA3';
                s(i).timeband = [1.231];  % 
                s(i).freqband = [7.812];
                s(i).label = 'Animal L Delay Cat Alpha2';
                
            i=i+1;
                s(i).label_short = 'LDA4';
                s(i).timeband = [.721];  % 
                s(i).freqband = [1.953];
                s(i).label = 'Animal L Delay Cat Delta';

            % % % % % % % % % % % % % % % % % % % % Scheme B % % % % % % % % % % % % % % % % % % % % 
            i=i+1;
%                 s(i).label_short = 'LSB1';
%                 s(i).timeband = [.541];  % 
%                 s(i).freqband = [21.48];
%                 s(i).label = 'Animal L Sample Goc Beta';
                
                s(i).label_short = 'LSB1';
                s(i).timeband = [.401];  % 
                s(i).freqband = [22.46];
                s(i).label = 'Animal L Sample Goc Beta';

        end
        if ~opts_exclude.excludeO
            % % % % % % % % % % % % % % % % % % % % Scheme A % % % % % % % % % % % % % % % % % % % % 
            i=i+1;
                s(i).label_short = 'OSA1';
                s(i).timeband = [.301];   % 
                s(i).freqband = [17.58];
                s(i).label = 'Animal O Sample Dog Beta';
                
            i=i+1; 
                s(i).label_short = 'OSA2';
                s(i).timeband = [.301];   % 
                s(i).freqband = [5.859];
                s(i).label = 'Animal O Sample Cat Theta';

            i=i+1; 
                s(i).label_short = 'ODA1';
                s(i).timeband = [.481];   % 
                s(i).freqband = [20.51];
                s(i).label = 'Animal O Delay Cat Beta';
                
            i=i+1; 
                s(i).label_short = 'ODA2';
                s(i).timeband = [.601];   % 
                s(i).freqband = [7.812];
                s(i).label = 'Animal O Delay Cat Alpha';
                
                
            i=i+1; 
                s(i).label_short = 'ODA3';
                s(i).timeband = [.961];   % 
                s(i).freqband = [15.62];
                s(i).label = 'Animal O Delay Cat Beta';
                
            i=i+1; 
                s(i).label_short = 'ODA4';
                s(i).timeband = [.871];   % 
                s(i).freqband = [3.906];
                s(i).label = 'Animal O Delay Dog Theta';

            % % % % % % % % % % % % % % % % % % % % Scheme B % % % % % % % % % % % % % % % % % % % % 
            
            i=i+1;
                s(i).label_short = 'OSB1';
                s(i).timeband = [.421];   % 
                s(i).freqband = [11.72];
                s(i).label = 'Animal O Sample Tad Alpha';
                
            i=i+1;
                s(i).label_short = 'OSB2';
                s(i).timeband = [.421];   % 
                s(i).freqband = [3.906];
                s(i).label = 'Animal O Sample Goc Theta';

            i=i+1;
                s(i).label_short = 'ODB1';
                s(i).timeband = [.901];   % 
                s(i).freqband = [5.859];
                s(i).label = 'Animal O Delay Tad Theta';

            i=i+1;
                s(i).label_short = 'ODB2';
                s(i).timeband = [1.291];   % 
                s(i).freqband = [17.58];
                s(i).label = 'Animal O Delay Tad Beta';

            i=i+1;
                s(i).label_short = 'ODB3';
                s(i).timeband = [1.441];   % 
                s(i).freqband = [11.72];
                s(i).label = 'Animal O Delay Goc Alpha';

            i=i+1;
                s(i).label_short = 'ODB4';
                s(i).timeband = [1.441];   % 
                s(i).freqband = [1.953];
                s(i).label = 'Animal O Delay Tad Delta';
                
            i=i+1;
                s(i).label_short = 'ODB5';
                s(i).timeband = [1.2];   % 
                s(i).freqband = [40];
                s(i).label = 'Animal O Delay Random Gamma';
                
            i=i+1;
                s(i).label_short = 'OSB3';
                s(i).timeband = [.541];   % 
                s(i).freqband = [18.55];
                s(i).label = 'Animal O Sample Goc and Tad (Both) Beta';
        end
    

    end


end