

%%% simu partially marker
%%% 1 2 3 3
%%% 1 2 3 4
%%% 1 2 2 2


%%%%产生前两个marker

tempr=3/4*(1-exp(-4*( tval.map(1) /100)/3));
tempphi=tempr^2/(9*(1-tempr)^2+tempr^2);
temppsi=tempr/(3-2*tempr);

Aeqt2(1,:)=[1,1,1,1,0,0,0,0,0];
Aeqt2(2,:)=[1,1,0,0,1,1,0,0,0];
Aeqt2(3,:)=[0,1,1/2,1,1/2,1,tempphi,(1+temppsi)/2,1];
Aeqt2(4,:)=1;
Beqt2=[tval.dr(1);tval.dr(2);tempr;1];

rg=rand(1,9);
rg=rg/sum(rg);
tempg=rand(1,9);

tempg=fmincon(@rand9g,tempg,[],[],Aeqt2,Beqt2,zeros(1,9),ones(1,9),[],optionsSD);
tempfreq=tempg(Mmtrx{1})./Mmtrx{2};
% SM.samp=zeros(10,10);
for is=1:samplesize
    tempj=rand(1);tempsum=0;
    for im1=1:10
        for im2=1:10
            tempsum=tempsum+tempfreq(im1,im2);
            if tempj<=tempsum
%                 SM.samp(im1,im2)=SM.samp(im1,im2)+1;
                SM.marker(is,1:2)=[im1,im2];
                tempj=9999;
            end
        end
    end
end

clear tempr tempphi temppsi tempg tempfreq is im1 im2 tempsum tempj





%%%% 2点模型产生其他相关marker

for i=3:(length(tval.map)+1)
    SM.marker(:,i)=SimuMarker2p(SM.marker(:,i-1),tval.dr(i-1:i),tval.map(i-1));
    if i==tval.loc
        SM.pai=SM_pai.pai;
        SM.paifreq=SM_pai.freq;
    end
end

%%%% 改变maker 3,5 phase
% 3 1234---1233
% 5 1234---1222
phase3{1}=1;phase3{2}=2;phase3{3}=[3,4,10];phase3{4}=5;phase3{5}=[6,7];phase3{6}=[8,9];
phase5{1}=1;phase5{2}=[2,3,4,8,9,10];phase5{3}=[5,6,7];

for ip1=1:length(phase3)
    for ip2=1:length(phase3{ip1})
        SM.marker(SM.marker(:,3)==phase3{ip1}(ip2),3)=ip1;
    end
end

for ip1=1:length(phase5)
    for ip2=1:length(phase5{ip1})
        SM.marker(SM.marker(:,5)==phase5{ip1}(ip2),5)=ip1;
    end
end

clear ip1 ip2

%%%% 产生表型

SM.phen=normrnd(tval.u(SM.marker(:,tval.loc)),sqrt(tval.sigma2))';














