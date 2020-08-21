
%%%d

dd=d([1,17:19,4:13,15,21:30],:);
dmean=mean(dd,2);
for i=1:length(dd(:,1))
    dstd(i)=std(dd(i,:));
end
dstd=dstd';


