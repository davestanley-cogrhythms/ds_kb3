

function normalize_fonts

 h= gcf; %get handle to current figure 
 axes = findobj(h.Children,'type','axes'); %get axes in figure

 for ax = axes(:)' %loop over axes
   ax.FontUnits = 'normalized'; %set units to normalized
 end

end