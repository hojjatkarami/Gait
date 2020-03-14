function Click_recon(cb,evendata,P,i_sub,i_trial,i_mus,SYN)
%cb is the handle of the axes that was clicked
%click on the whitespace within and axes and not on the line object
%copy the axes object to the new figure
% set(cb,'XTick',1:length(label))

% ClickedSubPlot = copyobj(cb,figure);
%for the new figure assign the ButtonDownFcn to empty
% set(ClickedSubPlot,'ButtonDownFcn',[]);
figure;
k=1;
for i_syn = SYN-2:SYN+2
    subplot(5,2,2*k-1); hold on
    plot(P(i_sub).EMG.Right(i_trial).Mn(i_mus,:))
    gof = P(i_sub).Synergy.EMG.Right(i_trial).gof;
    plot(P(i_sub).Synergy.EMG.Right(i_trial).syn(i_syn).M_rec(i_mus,:))
    title(['syn: ',num2str(i_syn),' vaf: ',num2str(gof.vaf(i_syn,i_mus),'%.2f')])
    k=k+1;
end
%resize the axis to fill the figure
% set(ClickedSubPlot, 'Position', get(0, 'DefaultAxesPosition'));
% h=get(ClickedSubPlot,'Children');
% xticks(1:length(label))
% xtickangle(45)
% xticklabels(label)
