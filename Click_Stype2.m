function Click_Stype2(cb,evendata)
gaitPhase = [12 31 50 62 75 87];
global idxcorr1 idxcorr2 idxcorr2m idxcorr3 xcorrData2 subplotID axesHandle currentFigure_S
f=gcf;
clickedFigNum=f.Number;
mouseButton = evendata.Button;
ClickedID = str2num(cb.Tag);
ClickedSubPlot = copyobj(cb,figure(222));
%for the new figure assign the ButtonDownFcn to empty
set(ClickedSubPlot,'ButtonDownFcn',[]);
%resize the axis to fill the figure
set(ClickedSubPlot, 'Position', get(0, 'DefaultAxesPosition'));
h=get(ClickedSubPlot,'Children');
switch mouseButton
    case 1  % compare two activation profiles
%         delete(figure(222))
    case 2  % compare emg kin and emgkin profiles
        if idxcorr2==idxcorr2m    
            figure(1000)
            clf(1000)
            idxcorr2=0;
        end
        if idxcorr2==0
%             answer = inputdlg('Number of Inputs?','Input',1,{'3'});
%             idxcorr2m = str2num(cell2mat(answer));
            idxcorr2m=2;
        end
        idxcorr2 = idxcorr2+1;
        axesHandle{idxcorr2} = h;
        xcorrData2{idxcorr2} = get(h(end), 'YData');
            figure(1000)
            axS = subplot(idxcorr2m,2,4);  copyobj(h,axS)
            axW = subplot(idxcorr2m,2,2*idxcorr2-1);
            
            figW=figure(clickedFigNum-1);
            allAxes = figW.Children;
            selectedAxes = allAxes(length(allAxes)+1-ClickedID);
            figure(clickedFigNum);
            figure(1000);
            copyobj(selectedAxes.Children,axW);
            xticks(axW,selectedAxes.XTick)
            xticklabels(axW,selectedAxes.XTickLabel)
            xtickangle(axW,45)
            if idxcorr2~=1
                for i=1:length(gaitPhase)
                   line(axS,[1 1]*gaitPhase(i),[0 8],'color','black');                     
                end
                [s1 s2, t21,befCorr,aftCorr] = xcorr0602(xcorrData2{idxcorr2-1}',xcorrData2{idxcorr2}');
                subplot(idxcorr2m,2,2*idxcorr2); hold on;
                plot(s2,'k--')
                title(['lag:',num2str(t21),' /Corr: ',num2str(befCorr),' - ',num2str(aftCorr)])

            end
        
        delete(figure(222));
    case 3  % just open it
        
end






