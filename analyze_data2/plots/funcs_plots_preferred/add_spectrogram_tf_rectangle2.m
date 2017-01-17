

function add_spectrogram_tf_rectangle2 (center_T,center_F,fullband_T,fullband_F,mylabel,mycolor, mylinestyle, show_text, mylegend)

    if nargin < 8
        show_text = 1;
    end
    
    if nargin < 7
        mylinestyle = '-';
    end
    
    if nargin < 6
        mycolor = 'w';
    end
    
    if nargin < 5
        mylabel = '';
    end
    
    tstart = [center_T - fullband_T/2];
    fstart = [center_F - fullband_F/2];
    hold on;
    rectangle('Position',[tstart,fstart,fullband_T,fullband_F],'EdgeColor',mycolor,'LineWidth',4,'LineStyle',mylinestyle);
    plot([tstart + fullband_T/2], [fstart + fullband_F/2],[mycolor 'x'],'MarkerSize',10,'LineWidth',2);
    
    if exist('mylegend','var')
        showtext_flag = any(strcmp_substr(mylegend,{'Cat','Dog'})) && any(strcmp_substr(mylabel,{'Cat','Dog'})) || ...
               any(strcmp_substr(mylegend,{'Goc','Tad'})) && any(strcmp_substr(mylabel,{'Goc','Tad'}));
    else
        showtext_flag = 1;
    end
    
    if show_text && showtext_flag
            text(tstart+.75*fullband_T,fstart+fullband_F,mylabel,'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',12,'Color','k');
    end            
end


% function b = strcontains(mylabel,mytext)
%     if iscell(mytext)
%         for i = 1:length(mytext)
%             b = b || ~isempty(strfind(mylabel{i},mytext));      % If any Any are true...
%         end
%     else
%         b = ~isempty(strfind(mylabel,mytext));
%     end
% end