function [s1, s2, t21,befCorr,aftCorr] = xcorr0602(s1,s2)                 
        
                   %
            [C21,lag21] = xcorr(s2,s1,50);
            C21 = C21/max(C21);
            %
            [M21,I21] = max(C21);
            t21 = lag21(I21);

            %
%             figure(111)
%             ax(2) = subplot(3,1,3);
%             plot(lag21,C21,[t21 t21],[-0.5 1],'r:')
%             text(t21+100,0.5,['Lag: ' int2str(t21)])
%             ylabel('C_{21}')
%             % axis tight
%             title(['Cross-Correlations, Lag: ',int2str(t21)])
% 
%             xlabel('Samples')
            %
            befCorr = num2str(corr(s1,s2),'%.2f');
            if t21>0
%                 s1 = delayseq(s1,t21,1);
%                 subplot(3,1,1)
%                 hold on
%                 plot(s1,'k--')
                
                s1=s1(1:end-t21);
                s2=s2(t21+1:end);
%                 s1=s1(1:end-t21);
%                 s2=s2(t21+1:end);
            elseif t21<0
                s2=s2(1:end+t21)
                s1=s1(-t21+1:end)
%                 s2 = delayseq(s2,-t21,1);
%                 subplot(3,1,2)
%                 hold on
%                 plot(s2,'k--')
%                 s2=s2(1:end-(-t21));
%                 s1=s1(-t21+1:end);

            end
            aftCorr = num2str(corr(s1,s2),'%.2f');
%             suptitle([befCorr,' -- ',aftCorr]) 


end