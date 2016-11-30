

function plot_venn(varargin)
% plot_venn(A,B,mylabels,work_in_percents,fixed_sizes)   - 3-way venn
% plot_venn(A,B,C,mylabels,work_in_percents,fixed_sizes) - 2-way venn

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

if nargin == 2
    plot_venn2(varargin{:});
end

if nargin >= 3
    if iscell(varargin{3})
        plot_venn2(varargin{:});
    else
        plot_venn3(varargin{:});
    end
end

end

function plot_venn3(A,B,C,mylabels,work_in_percents,fixed_sizes)

    if nargin < 5; work_in_percents = 0;
    end
    
    if nargin < 6; fixed_sizes = 0;
    end
    
    if work_in_percents; scaling_factor = sum(A | B | C) / 100;
        Nfiles = 100;
    else scaling_factor = 1;
        Nfiles = length(A);
    end
    
    pA = sum(A & ~B & ~C)/scaling_factor;     %A = clipping
    pB = sum(~A & B & ~C)/scaling_factor;     %B = 60 Hz
    pC = sum(~A & ~B & C)/scaling_factor;     %C = lowfiring
    
    
    pAB = sum(A & B & ~C)/scaling_factor;
    pAC = sum(A & ~B & C)/scaling_factor;
    pBC = sum(~A & B & C)/scaling_factor;
    
    pABC = sum(A & C & B)/scaling_factor;
    
    
    if nargin < 4
        venn([pA,pB,pC],[pAB,pAC,pBC,pABC]);
    else
        venn_labeled(mylabels,work_in_percents,fixed_sizes,[pA,pB,pC],[pAB,pAC,pBC,pABC]);
    end

end

function plot_venn2(A,B,mylabels,work_in_percents,fixed_sizes)

    if nargin < 4; work_in_percents = 0;
    end
    if nargin < 5; fixed_sizes = 0;
    end
    
    if work_in_percents; scaling_factor = sum(A | B ) / 100;
        Nfiles = 100;
    else scaling_factor = 1;
        Nfiles = length(A);
    end
    
    pA = sum(A & ~B)/scaling_factor; %A = clipping
    pB = sum(B & ~A)/scaling_factor;    %B = 60 Hz
    
    
    pAB = sum(A & B)/scaling_factor;
    
    
    if nargin < 3
        venn([pA,pB],[pAB]);
    else
        venn_labeled(mylabels,work_in_percents,fixed_sizes,[pA,pB],[pAB]);
    end
end


