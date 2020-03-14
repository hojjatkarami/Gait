function Click_showXtick(cb,evendata,label)

if isempty(cb.XTick)
    set(cb,'XTick',1:length(label))
    set(cb,'XTickLabelRotation',45)
    set(cb,'xticklabels',label)
else
    set(cb,'XTick',[])

    
end