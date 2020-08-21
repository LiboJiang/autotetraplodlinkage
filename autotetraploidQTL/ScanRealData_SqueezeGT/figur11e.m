

clear
load SQdata

load ScanResult
lr=result(2,:);
plot(lr,'k');

figure(gcf)

hold on;


x=1:300;
y=387*ones(1,300);
plot(x,y,'k-.');
hold on;


for i=1:length(scan.chr)                        %%mk---------------0~413.7
    x=0:0.1:20;
    y=scan.chr(i)*2*ones(1,length(x));
    plot(y,x,'k','linewidth',1.5);
    hold on;
end


group=[47.5,85,106.5];

for i=1:length(group)                        %%mk---------------0~413.7
    x=0:0.1:1000;
    y=group(i)*2*ones(1,length(x));
    plot(y,x,'k-.','linewidth',1);
    hold on;
end

QTL=263;



x=0:0.1:(fix(result(2,QTL)));
y=QTL*ones(1,length(x));
plot(y,x,'k','linewidth',1);
plot(QTL,25,'kV');

ylabel('Log-liklihood Ratio (LR)');
hold on
plot(0,1000)
hold on
axis tight






