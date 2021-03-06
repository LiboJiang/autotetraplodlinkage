
%%%
%%% mapping QTL simulation
%%%
%%%%%  sum(eval.pai(1:4))
%%%%%  sum(eval.pai([1,2,5,6]))


clc;
clear;

warning off all

global eval Mmtrx rg samp MQMmtrx

optionsEP=optimset('Algorithm','interior-point','DiffMinChange',1e-12,'MaxIter',50,'MaxFunEvals',7000,'Display','off');
optionsSD=optimset('Algorithm','interior-point', 'Display','off');%,'MaxFunEvals',5000);

rand('state',sum(100*clock));

%%%% Set parameter
%%%%%%%
tval.h2=0.2;
samplesize=200;
%%%%%%%
tval.ovau=1;
tval.add=[0.4,0.5,0.6];
tval.dom=[0.3,0.4,0.5,0.4,0.5,0.5];      %12,13,14,23,24,34

tval.loc=4;
tval.map=[10,10,12,8,10,10,10];
tval.dr=[0.2,0.2,0.2,0.2,0.2,0.2,0.2,0.2];

simutime=500;
savename=strcat('SimuResult',num2str(samplesize),'_',num2str(tval.h2*100));

loaddataSM;

SetScan.start=14;
SetScan.end=50;
SetScan.dens=2;
SimuResult=[];

fprintf('Simulation Strat  Samplesize=%d   H2=%0.2f  Scan:%d---%d\n\n',samplesize,tval.h2,SetScan.start,SetScan.end)
  
simu=0;

while simu<=simutime
err=0;
simu=simu+1;
result=[];
simudata;

scan.point=SetScan.start;
scan.mk=0;scan.sideL=0;scan.sideR=-1;

  while 1
    while scan.point>scan.sideR&&scan.sideR~=0
        scan.mk=scan.mk+1;scan.sideL=scan.chr(scan.mk);scan.sideR=scan.chr(scan.mk+1);scan.resamp=1;
    end
    if scan.sideR==0 
         break;
    end
   
    eval.r=[3/4*(1-exp(-4*((scan.point-scan.sideL)/100)/3)),3/4*(1-exp(-4*((scan.sideR-scan.point)/100)/3))];
    eval.r(3)=eval.r(1)+eval.r(2)-4/3*eval.r(1)*eval.r(2);
    eval.phi(1)=eval.r(1)^2/(9-18*eval.r(1)+10*eval.r(1)^2);
    eval.phi(2)=eval.r(2)^2/(9-18*eval.r(2)+10*eval.r(2)^2);
    eval.psi(1)=eval.r(1)/(3-2*eval.r(1));
    eval.psi(2)=eval.r(2)/(3-2*eval.r(2));
    
    if scan.resamp
        for i=1:10
            for j=1:10
                samp.PHmtrx{i,j}=[];
            end
        end
        samp.Nmtrx=zeros(10,10);
        for i=1:length(SM.marker(:,1)) 
            samp.Nmtrx(SM.marker(i,scan.mk),SM.marker(i,scan.mk+1))=samp.Nmtrx(SM.marker(i,scan.mk),SM.marker(i,scan.mk+1))+1;
            samp.PHmtrx{SM.marker(i,scan.mk),SM.marker(i,scan.mk+1)}=[samp.PHmtrx{SM.marker(i,scan.mk),SM.marker(i,scan.mk+1)},SM.phen(i)];
        end
        eval.pai(1)=(samp.Nmtrx(1,1)+samp.Nmtrx(2,2)+samp.Nmtrx(3,3)+samp.Nmtrx(4,4))/samplesize;
        eval.pai(2)=(samp.Nmtrx(1,2)+samp.Nmtrx(1,3)+samp.Nmtrx(1,4)+samp.Nmtrx(2,1)+samp.Nmtrx(2,3)+samp.Nmtrx(2,4)+samp.Nmtrx(3,1)+samp.Nmtrx(3,2)+samp.Nmtrx(3,4)+samp.Nmtrx(4,1)+samp.Nmtrx(4,2)+samp.Nmtrx(4,3))/samplesize;
        eval.pai(3)=(samp.Nmtrx(1,5)+samp.Nmtrx(1,6)+samp.Nmtrx(1,7)+samp.Nmtrx(2,5)+samp.Nmtrx(2,8)+samp.Nmtrx(2,9)+samp.Nmtrx(3,6)+samp.Nmtrx(3,8)+samp.Nmtrx(3,10)+samp.Nmtrx(4,7)+samp.Nmtrx(4,9)+samp.Nmtrx(4,10))/samplesize;
        eval.pai(4)=(samp.Nmtrx(1,8)+samp.Nmtrx(1,9)+samp.Nmtrx(1,10)+samp.Nmtrx(2,6)+samp.Nmtrx(2,7)+samp.Nmtrx(2,10)+samp.Nmtrx(3,5)+samp.Nmtrx(3,7)+samp.Nmtrx(3,9)+samp.Nmtrx(4,5)+samp.Nmtrx(4,6)+samp.Nmtrx(4,8))/samplesize;
        eval.pai(5)=(samp.Nmtrx(5,1)+samp.Nmtrx(5,2)+samp.Nmtrx(6,1)+samp.Nmtrx(6,3)+samp.Nmtrx(7,1)+samp.Nmtrx(7,4)+samp.Nmtrx(8,2)+samp.Nmtrx(8,3)+samp.Nmtrx(9,2)+samp.Nmtrx(9,4)+samp.Nmtrx(10,3)+samp.Nmtrx(10,4))/samplesize;
        eval.pai(6)=(samp.Nmtrx(5,3)+samp.Nmtrx(5,4)+samp.Nmtrx(6,2)+samp.Nmtrx(6,4)+samp.Nmtrx(7,2)+samp.Nmtrx(7,3)+samp.Nmtrx(8,1)+samp.Nmtrx(8,4)+samp.Nmtrx(9,1)+samp.Nmtrx(9,3)+samp.Nmtrx(10,1)+samp.Nmtrx(10,2))/samplesize;
        eval.pai(7)=(samp.Nmtrx(5,5)+samp.Nmtrx(6,6)+samp.Nmtrx(7,7)+samp.Nmtrx(8,8)+samp.Nmtrx(9,9)+samp.Nmtrx(10,10))/samplesize;
        eval.pai(8)=(samp.Nmtrx(5,6)+samp.Nmtrx(5,7)+samp.Nmtrx(5,8)+samp.Nmtrx(5,9)+samp.Nmtrx(6,5)+samp.Nmtrx(6,7)+samp.Nmtrx(6,8)+samp.Nmtrx(6,10)+samp.Nmtrx(7,5)+samp.Nmtrx(7,6)+samp.Nmtrx(7,9)+samp.Nmtrx(7,10)...
                +samp.Nmtrx(8,5)+samp.Nmtrx(8,6)+samp.Nmtrx(8,9)+samp.Nmtrx(8,10)+samp.Nmtrx(9,5)+samp.Nmtrx(9,7)+samp.Nmtrx(9,8)+samp.Nmtrx(9,10)+samp.Nmtrx(10,6)+samp.Nmtrx(10,7)+samp.Nmtrx(10,8)+samp.Nmtrx(10,9))/samplesize;
        eval.pai(9)=(samp.Nmtrx(5,10)+samp.Nmtrx(6,9)+samp.Nmtrx(7,8)+samp.Nmtrx(8,7)+samp.Nmtrx(9,6)+samp.Nmtrx(10,5))/samplesize;    
        
        scan.resamp=0;
    end
    
    
    Aeq=zeros(13,59);
    for i=1:9
        Aeq(i,coefx.pai{i})=1;     %pai
    end
    Aeq(10,coefx.r{1,1})=1;Aeq(10,coefx.r{1,2})=1/2;Aeq(10,coefx.r{1,3})=eval.phi(1);Aeq(10,coefx.r{1,4})=1/2*(1+eval.psi(1));  %r1
    Aeq(11,coefx.r{2,1})=1;Aeq(11,coefx.r{2,2})=1/2;Aeq(11,coefx.r{2,3})=eval.phi(2);Aeq(11,coefx.r{2,4})=1/2*(1+eval.psi(2));  %r2
    Aeq(12,[5:11,19:25,33:39,52,55:59])=1;Aeq(12,[12:18,26:32,40,41,44,49,50,51])=1/2;Aeq(12,[45,46,47,48,53,54])=3/4;   %r3
    Aeq(12,42)=Aeq(12,42)+eval.phi(1);Aeq(12,49)=Aeq(12,49)+eval.phi(1)/2;
    Aeq(12,42)=Aeq(12,42)+eval.phi(2);Aeq(12,50)=Aeq(12,50)+eval.phi(2)/2;
    Aeq(12,42)=Aeq(12,42)-2*eval.phi(1)*eval.phi(2);
    Aeq(12,43)=Aeq(12,43)+eval.psi(1);Aeq(12,[50,51])=Aeq(12,[50,51])+eval.psi(1)/2;Aeq(12,52)=Aeq(12,52)-eval.psi(1)/2;Aeq(12,50)=Aeq(12,50)-eval.phi(2)*eval.psi(1);            
    Aeq(12,43)=Aeq(12,43)+eval.psi(2);Aeq(12,[49,51])=Aeq(12,[49,51])+eval.psi(2)/2;Aeq(12,52)=Aeq(12,52)-eval.psi(2)/2;Aeq(12,49)=Aeq(12,49)-eval.phi(1)*eval.psi(2);
    Aeq(12,43)=Aeq(12,43)-2*eval.psi(1)*eval.psi(2);Aeq(12,51)=Aeq(12,51)-eval.psi(1)*eval.psi(2);Aeq(12,52)=Aeq(12,52)+eval.psi(1)*eval.psi(2);
    Aeq(13,:)=1;
    Beq=[eval.pai';eval.r';1];
    
    %%%% est para
    
    eval.u=ones(1,10);eval.sigma2=3;eval.g=(ones(1,59)/59)';
    Ou=zeros(1,10);Osigma2=0;
    CRtime=0;


    while sum(abs([eval.u,eval.sigma2]-[Ou,Osigma2])>0.001)>0
        [eval.g,fval,exitflag,output]=fmincon(@ComputeL1,eval.g,[],[],Aeq,Beq,zeros(1,59),ones(1,59),[],optionsEP);
        Ou=eval.u;Osigma2=eval.sigma2;
        part_Em;
        CRtime=CRtime+1;
        fprintf('Simu/Scan %d/%.1f  CRtime=%d--U:(%6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f)    Sigma2=%6.3f   #emtime=%d\n',simu,scan.point,CRtime,eval.u,eval.sigma2,emtime)
    end
    
    L1=-ComputeL1(eval.g);
    L0=ComputeL0(mean(SM.phen),var(SM.phen));
    LR=-2*(L0-L1);
        
    eval.dr(1)=sum(eval.pai([1,2,3,4]));
    eval.dr(2)=sum(eval.pai([1,2,5,6]));
    eval.dr(3)=sum(eval.g([1,2,5,6,7,12,13,14,19,20,21,26,27,28,33,34,35,40,41,45,46,47,48,55,56]));
    eval.par(1)=sum(eval.u(1:4))/4;                                             % over all u
    eval.par(2:4)=eval.u(1:3)-eval.par(1);                                      % add effect a1 a2 a3
    eval.par(5)=eval.u(5)-sum(eval.par([1,2,3]));                               % d12
    eval.par(6)=eval.u(6)-sum(eval.par([1,2,4]));                               % d13
    eval.par(7)=eval.u(7)-eval.par(1)+eval.par(3)+eval.par(4);                  % d14
    eval.par(8)=eval.u(8)-sum(eval.par([1,3,4]));                               % d23
    eval.par(9)=eval.u(9)-eval.par(1)+eval.par(2)+eval.par(4);                  % d24
    eval.par(10)=eval.u(10)-eval.par(1)+eval.par(2)+eval.par(3);                % d34
    result=[result,[scan.point;LR;0;eval.u';0;eval.sigma2;0;eval.dr';0;eval.par']];
    fprintf('**simu/scan %d/%.1f  LR=%.2f   U:(%6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f)    Sigma2=%6.3f\n',simu,scan.point,LR,eval.u,eval.sigma2) 
    
    
    scan.point=scan.point+SetScan.dens;
    if scan.point>SetScan.end
        break;
    end
    
    clear CRtime L1 L0 LR Osigma2 Ou ou osigma2 etime exitflag fval i j l tempsg tempsum tempu1 tempu2 output freqP
    
  end
  
  
  result1=result(:,result(2,:)==max(result(2,:)));
  fprintf('\n*****simu %d has been done\n',simu);
  fprintf('simu/scan %d/%.1f    Max LR=%.2f  U:%6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f   Sigma2=%6.3f\n',simu,result1([1,2,4:13,15]));
  fprintf('Dr: %6.3f  %6.3f  %6.3f\n',result1(17:19));
  fprintf(' u: %6.3f\n a: %6.3f %6.3f %6.3f \n d: %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f\n\n',result1(21:30));

   
  SimuResult(:,simu)=result1;
  save(savename,'SimuResult','samplesize','tval');
  clear result1 result SM Aeqt2 Beqt2
  
end

    

     








