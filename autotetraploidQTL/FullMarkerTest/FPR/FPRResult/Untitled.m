
n=3;

permu=FPRpermu(n,:);
fpr1=fpr(n,:);

permu=sort(permu);
p1=permu(length(permu)-4);   %---0.005
p2=permu(length(permu)-9);   %---0.01
p3=permu(length(permu)-49);  %---0.05

sum((fpr1-p1)>0)
sum((fpr1-p2)>0)
sum((fpr1-p3)>0)