function [Name, partitionLine] = muscleName(No)
Name={};
    allNames = {'TA' ,'PL', 'SOL', 'GC' 'GC_M' ...
                'RF','MH' 'HM','LH', 'VL' ,'VM','GMED' 'G_M' ,'GMAX' 'GLUT_MAX'...
                'IP' 'IPS' ,'AD' ,'TFL'...
                'IC' ,'LG' 'RG','RA' 'RECT_AB', 'EO'};
    allNo = {'1' ,'2', '3', '4' '4' ...
             '5', '6' '6' ,'7', '8' ,'9','10' '10' ,'11' '11'...
             '12' '12', '13' ,'14'...
             '15' ,'16' '16','17' '17', '18'};
for i=1:length(No)
    id = find(ismember(allNo,num2str(No(i))));    
    Name{i} = allNames{id(1)};
end


