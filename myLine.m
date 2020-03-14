function myLine(type,x,y)

switch type
    case 'v'
        for i = 1:length(x)
           line([x(i) x(i)],y) 
        end
        
    case 'h'
        for i = 1:length(x)
           line(x,[y(i) y(i)]) 
        end
end

