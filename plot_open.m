function plot_open(cb,evendata)
%cb is the handle of the axes that was clicked
%click on the whitespace within and axes and not on the line object
ClickedSubPlot = copyobj(cb,figure);
%for the new figure assign the ButtonDownFcn to empty
set(ClickedSubPlot,'ButtonDownFcn',[]);

%resize the axis to fill the figure
set(ClickedSubPlot, 'Position', get(0, 'DefaultAxesPosition'));
h=get(ClickedSubPlot,'Children');
    

