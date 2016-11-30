



function [h,hl] = plott_handles(fname,varargin)
%     [h,hl] = plott_handles(fname,fsubplot,varargin)
%     Takes in a cell array of function handles and plots them in a bunch of
%     subplots on a single figure
%     FORMS
%         plott_handles(fname,varargin)
%         plott_handles(fname,fsubplot,varargin)
%     
%     INPUTS
%         fname - cell array of function handles to plot
%         fsubplot - subplot function to use. Default is @subplotsq
%         varargin - these inputs are just forwarded on to each of the fnames
%     OUTPUTS
%         h = collection of subplot object handles
%         hl = line handles from each subplot
% 
%     Example 1
%     load audio48
%     x = signal48kHz(1:round(end/10));
%     fnames = {@plot,@loglog,@plott_psd,@(t,x) plott_spect(t,x,'mode',1,'logplot',1)};
%     t = 1:length(x); t=t/Fs48;
%     figure; plott_handles(fnames,t,x);
%     figure; plott_handles(fnames,@subplotrows,t,x); % Plot with different subplot arrangement
% 
%     CONTACT: David Stanley, Boston University (stanleyd@bu.edu, https://github.com/davestanley)
% 


    % Import into workspace
    if isa(varargin{1},'function_handle')
        fsubplot = varargin{1};
        args = varargin(2:end);
    else
        fsubplot=@subplotsq;
        args = varargin;
    end
    
    if isempty(fsubplot)
        fsubplot=@subplotsq;
    end
    
    
    h=[];
    if iscell(fname)
        Nplots = length(fname);
        for j=1:Nplots
            h(j) = fsubplot(Nplots,j);
            hl{j} = fname{j}(args{:}); hold on; 
        end
    else
        hold on; hl = fname(args{:});
    end
    
end

