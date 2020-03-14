function [lb,ub] = myBootStrap(n,percent,func,M,rec)


eval(['[bootstat,~] = bootstrp(n,',func,',M,rec,0);']);
r2_sorted = sort(bootstat);
alpha = 1- percent/100;
lb = r2_sorted(round(alpha*n));
ub = r2_sorted(n-round(alpha*n));
