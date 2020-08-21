

%%% simu partially marker
%%% 1 0 0 0
%%% 1 2 2 2
%%% 1 0 0 0


%%%%产生前两个marker

    tval.r=[3/4*(1-exp(-4*((tval.map(1))/100)/3)),3/4*(1-exp(-4*((tval.map(2))/100)/3))];
    tval.r(3)=tval.r(1)+tval.r(2)-4/3*tval.r(1)*tval.r(2);
    tval.phi(1)=tval.r(1)^2/(9-18*tval.r(1)+10*tval.r(1)^2);
    tval.phi(2)=tval.r(2)^2/(9-18*tval.r(2)+10*tval.r(2)^2);
    tval.psi(1)=tval.r(1)/(3-2*tval.r(1));
    tval.psi(2)=tval.r(2)/(3-2*tval.r(2));

    Aeqs=zeros(4,59);

    Aeqs(1,coefx.r{1,1})=1;Aeqs(1,coefx.r{1,2})=1/2;Aeqs(1,coefx.r{1,3})=tval.phi(1);Aeqs(1,coefx.r{1,4})=1/2*(1+tval.psi(1));  %r1
    Aeqs(2,coefx.r{2,1})=1;Aeqs(2,coefx.r{2,2})=1/2;Aeqs(2,coefx.r{2,3})=tval.phi(2);Aeqs(2,coefx.r{2,4})=1/2*(1+tval.psi(2));  %r2
    Aeqs(3,[5:11,19:25,33:39,52,55:59])=1;Aeqs(3,[12:18,26:32,40,41,44,49,50,51])=1/2;Aeqs(3,[45,46,47,48,53,54])=3/4;   %r3
    Aeqs(3,42)=Aeqs(3,42)+tval.phi(1);Aeqs(3,49)=Aeqs(3,49)+tval.phi(1)/2;
    Aeqs(3,42)=Aeqs(3,42)+tval.phi(2);Aeqs(3,50)=Aeqs(3,50)+tval.phi(2)/2;
    Aeqs(3,42)=Aeqs(3,42)-2*tval.phi(1)*tval.phi(2);
    Aeqs(3,43)=Aeqs(3,43)+tval.psi(1);Aeqs(3,[50,51])=Aeqs(3,[50,51])+tval.psi(1)/2;Aeqs(3,52)=Aeqs(3,52)-tval.psi(1)/2;Aeqs(3,50)=Aeqs(3,50)-tval.phi(2)*tval.psi(1);            
    Aeqs(3,43)=Aeqs(3,43)+tval.psi(2);Aeqs(3,[49,51])=Aeqs(3,[49,51])+tval.psi(2)/2;Aeqs(3,52)=Aeqs(3,52)-tval.psi(2)/2;Aeqs(3,49)=Aeqs(3,49)-tval.phi(1)*tval.psi(2);
    Aeqs(3,43)=Aeqs(3,43)-2*tval.psi(1)*tval.psi(2);Aeqs(3,51)=Aeqs(3,51)-tval.psi(1)*tval.psi(2);Aeqs(3,52)=Aeqs(3,52)+tval.psi(1)*tval.psi(2);
    Aeqs(4,:)=1;
    Beqs=[tval.r';1];%;tval.pai];
    
%     for i=1:9
%         Aeqs(i,coefx.pai{i})=1;     %pai
%     end

  
    rg=rand(1,59);
    rg=rg/sum(rg);
    tempg=rand(1,59);

    tempg=fmincon(@rand59g,tempg,[],[],Aeqs,Beqs,zeros(1,59),ones(1,59),[],optionsSD);



% tempfreqs=tempg(Mmtrx{1})./Mmtrx{2};
% % SM.samp=zeros(10,10);
% for is=1:samplesize
%     tempj=rand(1);tempsum=0;
%     for im1=1:10
%         for im2=1:10
%             tempsum=tempsum+tempfreqs(im1,im2);
%             if tempj<=tempsum
% %                 SM.samp(im1,im2)=SM.samp(im1,im2)+1;
%                 SM.marker(is,1:2)=[im1,im2];
%                 tempj=9999;
%             end
%         end
%     end
% end
% 
% clear tempr tempphi temppsi tempg tempfreqs is im1 im2 tempsum tempj
% 
% 
% 
% 
% 
% %%%% 2点模型产生其他相关marker
% 
% for i=3:(length(tval.map)+1)
%     SM.marker(:,i)=SimuMarker2p(SM.marker(:,i-1),tval.dr(i-1:i),tval.map(i-1));
%     if i==tval.loc
%         SM.pai=SM_pai.pai;
%         SM.paifreqs=SM_pai.freqs;
%     end
% end
% 
% %%%% 改变maker 3,5 phase
% % 3 1234---1233
% % 5 1234---1222
% phase3{1}=1;phase3{2}=2;phase3{3}=[3,4,10];phase3{4}=5;phase3{5}=[6,7];phase3{6}=[8,9];
% phase5{1}=1;phase5{2}=[2,3,4,8,9,10];phase5{3}=[5,6,7];
% 
% for ip1=1:length(phase3)
%     for ip2=1:length(phase3{ip1})
%         SM.marker(SM.marker(:,3)==phase3{ip1}(ip2),3)=ip1;
%     end
% end
% 
% for ip1=1:length(phase5)
%     for ip2=1:length(phase5{ip1})
%         SM.marker(SM.marker(:,5)==phase5{ip1}(ip2),5)=ip1;
%     end
% end
% 
% clear ip1 ip2
% 
% %%%% 产生表型
% 
% SM.phen=normrnd(tval.u(SM.marker(:,tval.loc)),sqrt(tval.sigma2))';
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
