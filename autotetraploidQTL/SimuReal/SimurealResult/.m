
%%%d

dd=d([1,11:13,5:7,9],:);
for i=1:length(dd(1,:))
    dd(5:7,i)=sort(dd(5:7,i));
end
dmean=mean(dd,2);
for i=1:length(dd(:,1))
    dstd(i)=std(dd(i,:));
end
dstd=dstd';


