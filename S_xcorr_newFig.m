function S_xcorr_newFig(cb,evendata)
%cb is the handle of the axes that was clicked
%click on the whitespace within and axes and not on the line object
%copy the axes object to the new figure
global idxcorr xcorrData subplotID
if idxcorr==2
    
    figure(111)
    clf(111)
end
idxcorr =mod(idxcorr,2)+1;

ClickedSubPlot = copyobj(cb,figure(222+idxcorr));
%for the new figure assign the ButtonDownFcn to empty
set(ClickedSubPlot,'ButtonDownFcn',[]);

%resize the axis to fill the figure
set(ClickedSubPlot, 'Position', get(0, 'DefaultAxesPosition'));
h=get(ClickedSubPlot,'Children');


% copyobj(cb,subplotID{idxcorr+1})
if length(h)==2
    h = h(2);
end
figure(111)
subplot(3,1,idxcorr);
xcorrData{idxcorr} = get(h, 'YData');
plot(xcorrData{idxcorr},'Color',h.Color)

if idxcorr==2
    delete(223)
    delete(224)
                s1=xcorrData{1}'; s2=xcorrData{2}';

            %
            [C21,lag21] = xcorr(s2,s1,50);
            C21 = C21/max(C21);
            %
            [M21,I21] = max(C21);
            t21 = lag21(I21);

            %
            figure(111)
            ax(2) = subplot(3,1,3);
            plot(lag21,C21,[t21 t21],[-0.5 1],'r:')
            text(t21+100,0.5,['Lag: ' int2str(t21)])
            ylabel('C_{21}')
            % axis tight
            title(['Cross-Correlations, Lag: ',int2str(t21)])

            xlabel('Samples')
            %
            befCorr = num2str(corr(s1,s2),'%.2f');
            if t21>0
                s1 = delayseq(s1,t21,1);
                subplot(3,1,1)
                hold on
                plot(s1,'k--')
                s1=s1(1:end-t21);
                s2=s2(t21+1:end);
            elseif t21<0
                s2 = delayseq(s2,-t21,1);
                subplot(3,1,2)
                hold on
                plot(s2,'k--')
                s2=s2(1:end-(-t21));
                s1=s1(-t21+1:end);

            end
            aftCorr = num2str(corr(s1,s2),'%.2f');
            suptitle([befCorr,' -- ',aftCorr]) 


            end

end