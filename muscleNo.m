function No = muscleNo(name)
name = name(2:end);     % separate L or R char
LmuscleNames = {'LTA','LSOL','LGC_M','LRF','LMH','LVL','LIPS','LGLUT_MAX',...
    'LG_M','LVM','LTFL','LIC','LLG','LRECT_AB'};%%

allNames = {'TA' ,'PL', 'SOL', 'GC' 'GC_M' ...
            'RF','MH' 'HM','LH', 'VL' ,'VM','GMED' 'G_M' ,'GMAX' 'GLUT_MAX'...
            'IP' 'IPS' ,'AD' ,'TFL'...
            'IC' ,'LG' 'RG','RA' 'RECT_AB', 'EO'};
allNo = {'1' ,'2', '3', '4' '4' ...
        '5', '6' '6' ,'7', '8' ,'9','10' '10' ,'11' '11'...
        '12' '12', '13' ,'14'...
        '15' ,'16' '16','17' '17', '18'};
        
id = find(ismember(allNames,name));
if isempty(id)
    No = 0;
else
    No = str2num(cell2mat(allNo(id)));
end




