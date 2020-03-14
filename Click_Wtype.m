function Click_Wtype(cb,evendata,label)
%cb is the handle of the axes that was clicked
%click on the whitespace within and axes and not on the line object
%copy the axes object to the new figure
set(cb,'XTick',1:length(label))
ClickedSubPlot = copyobj(cb,figure);
%for the new figure assign the ButtonDownFcn to empty
set(ClickedSubPlot,'ButtonDownFcn',[]);

%resize the axis to fill the figure
set(ClickedSubPlot, 'Position', get(0, 'DefaultAxesPosition'));
h=get(ClickedSubPlot,'Children');
xticks(1:length(label))
xtickangle(45)
xticklabels(label)
