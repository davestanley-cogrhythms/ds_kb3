

function plot_venn_badfiles(bad_lowfiring, bad_clip, bad_60, bad_any)

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



    work_in_percents = 1;
    
    if work_in_percents; scaling_factor = length(bad_clip) / 100;
        Nfiles = 100;
    else scaling_factor = 1;
        Nfiles = length(bad_clip);
    end
    
    pA = sum(bad_clip)/scaling_factor; %A = clipping
    pB = sum(bad_60)/scaling_factor;    %B = 60 Hz
    pC = sum(bad_lowfiring)/scaling_factor;  %C = lowfiring
    
    
    pAB = sum(bad_clip & bad_60)/scaling_factor;
    pAC = sum(bad_clip & bad_lowfiring)/scaling_factor;
    pBC = sum(bad_60 & bad_lowfiring)/scaling_factor;
    
    pABC = sum(bad_60 & bad_lowfiring & bad_60)/scaling_factor;
    
    
    subplot(121);
    if sum(bad_lowfiring) < 1
        venn_labeled({'Clip','60Hz'},[pA,pB],[pAB]);
    else
        %venn_labeled({'Clip','60Hz'},[pA,pB],[pAB]);
        venn_labeled({'Clip','60Hz'},[pA,pB],[pAB]);
        %venn_labeled({'Clip','60Hz','LowFiring'},[pA,pB,pC],[pAB,pAC,pBC,pABC]);

    end
    fprintf('Bad files clipping %d; Bad files 60 Hz %d; Bad files lowfiring %d. \n', pA, pB,pC);
    
    
    pAny_bad = sum(bad_any)/scaling_factor;
    pAll_Good = sum(~bad_any)/scaling_factor;
    
    subplot(122);
%     venn([pAny_bad,pAll_Good,Nfiles],[0,pAny_bad,pAll_Good,0]);
    venn_labeled({'Bad','Good'},[pAny_bad,pAll_Good],[0]);
%     myrad = (sqrt(pAny_bad)+sqrt(pAll_Good)).^2;    % Circle having radius of pAny_bad+pAll_good
%     venn([myrad,pAny_bad,pAll_Good,],[pAny_bad,pAll_Good,0,0])    % 3rd circle for all
    


end


