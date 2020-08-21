

%%% simu partially marker
%%% 1 0 0 0
%%% 1 2 2 2
%%% 1 0 0 0


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
    Aeqs(4,coefx.dr{1})=1;
    Aeqs(5,coefx.dr{2})=1;
    Aeqs(6,coefx.dr{3})=1;
    Aeqs(7,:)=1;
    
    Beqs=[tval.r';tval.dr';1];%;tval.pai];
    
%     for i=1:9
%         Aeqs(i,coefx.pai{i})=1;     %pai
%     end

  
    rg=rand(59,1);
    rg=rg/sum(rg);
    tempg=rand(59,1);
    tval.g=fmincon(@rand59g,tempg,[],[],Aeqs,Beqs,zeros(1,59),ones(1,59),[],optionsSD);
    
    tval.pai=zeros(9,1);
    for i=1:9
        tval.pai(i)=sum(tval.g(coefx.pai{i}));     %pai
    end
    
    
    %%%%%  现在要利用 MQM 矩阵 计算m123
    
    SM.samp=zeros(10,10,10);
    tval.freqg=tval.g(MQMmtrx{1})./MQMmtrx{2};
    for is=1:samplesize
        tempj=rand(1);tempsum=0;
        for im1=1:10
            for im2=1:10
                for im3=1:10
                    tempsum=tempsum+tval.freqg(im1,im2,im3);
                    if tempj<=tempsum
                        SM.marker_F(is,1:3)=[im1,im2,im3];
                        SM.samp(im1,im2,im3)=SM.samp(im1,im2,im3)+1;
                        tempj=9999;
                    end
                end
            end
        end
    end
    

clear tempg tempfreqs is im1 im2 im3 tempsum tempj i




%%%% 改变phase
% 1 1234---1000
% 2 1234---1222
% 3 1234---1000
phase1{1}=[2,3,4,8,9,10];phase1{2}=[1,5,6,7];
phase2{1}=1;phase2{2}=[5,6,7];phase2{3}=[2,3,4,8,9,10];

for ii=1:2:3
    for ip1=1:length(phase1)
        for ip2=1:length(phase1{ip1})
            SM.marker_P(SM.marker_F(:,ii)==phase1{ip1}(ip2),ii)=ip1;
        end
    end
end

for ip1=1:length(phase2)
    for ip2=1:length(phase2{ip1})
        SM.marker_P(SM.marker_F(:,2)==phase2{ip1}(ip2),2)=ip1;
    end
end


clear ip1 ip2 ii

%%% 产生表型

SM.phen=normrnd(tval.u(SM.marker_P(:,tval.loc)),sqrt(tval.sigma2))';
SM.marker=SM.marker_P(:,[1,3]);

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
