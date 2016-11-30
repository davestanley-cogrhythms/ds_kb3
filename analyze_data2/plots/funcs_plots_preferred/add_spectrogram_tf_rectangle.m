

function add_spectrogram_tf_rectangle (center_T,center_F,fullband_T,fullband_F,mylabel,mycolor, mylinestyle, show_text, ctgs)

    
    
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
    
    if exist('ctgs','var')
        showtext_flag = ((ctgs == 1 || ctgs == 3) && (~isempty(strfind(mylabel,'Cat')) || ~isempty(strfind(mylabel,'Dog')))) || ...
               ((ctgs == 2 || ctgs == 4) && (~isempty(strfind(mylabel,'Goc')) || ~isempty(strfind(mylabel,'Tad')))) ;
    else
        showtext_flag = 1;
    end
    
    if show_text && showtext_flag
            text(tstart+.75*fullband_T,fstart+fullband_F,mylabel,'HorizontalAlignment','right','VerticalAlignment','bottom','FontSize',12,'Color','k');

    end            
end