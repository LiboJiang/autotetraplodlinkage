


n=2;

permu=d(n,1:1000);




permu=sort(permu);
p1=permu(length(permu)-1000*0.005+1);   %---0.005
p2=permu(length(permu)-1000*0.01+1);    %---0.01
p3=permu(length(permu)-1000*0.05+1);    %---0.05

fpr1=f(n,1:100);
sum((fpr1-p1)>0)
sum((fpr1-p2)>0)
sum((fpr1-p3)>0)