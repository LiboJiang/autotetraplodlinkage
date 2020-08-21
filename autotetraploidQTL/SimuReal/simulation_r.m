
%%%
%%%    mapping QTL simulation
%%%

%%%   partially informative
%%%   1 0 0 0
%%%   1 2 2 2
%%%   1 0 0 0 


clc;
clear;
clear global;


warning off all

global rg eval samp mdm samplesize scan Mmtrx Lphase SM 

optionsEP=optimset('Algorithm','interior-point','DiffMinChange',1e-12,'MaxIter',30,'MaxFunEvals',5000,'Display','off');
optionsSD=optimset('Algorithm','interior-point', 'Display','off');%,'MaxFunEvals',5000);

rand('state',sum(100*clock));

%%%% Set parameter
%%%%%%%
tval.h2=0.05;
samplesize=200;
%%%%%%%
tval.loc=2;
tval.map=[8,12];
tval.dr=[0.2,0.2,0.2];


tval.u=[140,160,180];
tempup1=tval.dr(tval.loc)/4;tempup2=(1-tval.dr(tval.loc))/6;
tval.up=[tempup1,3*tempup2,3*tempup1+3*tempup2];

tval.uu=sum(tval.u.*tval.up);
tval.sigma2=sum(((tval.u-tval.uu).^2.*tval.up))/tval.h2;

simutime=500;
savename=strcat('SimuReal',num2str(samplesize),'_',num2str(tval.h2*100));
qphase{1}=1;
qphase{2}=[2,3,4,8,9,10];
qphase{3}=[5,6,7];

loaddataSM;

SetScan.start=8;
SetScan.end=50;
SetScan.dens=17;
SimuResult=[];

fprintf('SIMU Samplesize=%d   H2=%0.2f  Scan:%d---%d\n\n',samplesize,tval.h2,SetScan.start,SetScan.end)
fprintf('Phase: 1 0 0 0\n')
fprintf('       1 2 2 2\n')
fprintf('       1 0 0 0\n\n')
  
simu=0;
clear tempup1 tempup2
while simu<=simutime
simu=simu+1;
result=[];
simudata_real;


scan.point=SetScan.start;
scan.mk=0;scan.sideL=0;scan.sideR=-1;

 while 1
    while scan.point>scan.sideR && scan.sideR~=0
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
        
        samp.MQMmtrxG=[];
        samp.MQMmtrxC=[];
        mdm(1)=2;mdm(2)=2;

        for i=1:mdm(1)
            for j=1:mdm(1)
                samp.PhenMtrx{i,j}=[];
            end
        end
        samp.Nmtrx=zeros(mdm(1),mdm(2));
        for i=1:length(SM.marker(:,1)) 
            samp.Nmtrx(SM.marker(i,scan.mk),SM.marker(i,scan.mk+1))=samp.Nmtrx(SM.marker(i,scan.mk),SM.marker(i,scan.mk+1))+1;
            samp.PhenMtrx{SM.marker(i,scan.mk),SM.marker(i,scan.mk+1)}=[samp.PhenMtrx{SM.marker(i,scan.mk),SM.marker(i,scan.mk+1)},SM.phen(i)];
        end

            
            %确定连锁相,计算pai
            tempLP=ConfirmLP(1);
            eval.pai=tempLP.epai;
           
            Mmtrx_LP=zeros(mdm(1),mdm(2));
            for i1=1:mdm(1)
                for i2=1:mdm(2)
                    Mmtrx_LP(i1,i2)=sum(eval.pai(tempLP.Mpai{1}{i1,i2})./tempLP.Mpai{2}{i1,i2});
                end
            end
            
            samp.Pai_Mmtrx=Mmtrx_LP;
                       
            for im1=1:2
                for im2=1:3
                    for im3=1:2
                        tempG{im1,im2,im3}=MQMmtrx{1}(tempLP.mtc{1,im1},qphase{im2},tempLP.mtc{2,im3});
                        tempC{im1,im2,im3}=MQMmtrx{2}(tempLP.mtc{1,im1},qphase{im2},tempLP.mtc{2,im3});
                    end
                end
            end
        
            for im1=1:2
                for im2=1:3
                    for im3=1:2
                        samp.MQMmtrxG{im1,im2,im3}=[];samp.MQMmtrxC{im1,im2,im3}=[];
                        for i1=1:length(tempG{im1,im2,im3}(:,1,1))
                            for i2=1:length(tempG{im1,im2,im3}(1,:,1))
                                for i3=1:length(tempG{im1,im2,im3}(1,1,:))
                                    samp.MQMmtrxG{im1,im2,im3}=[samp.MQMmtrxG{im1,im2,im3},tempG{im1,im2,im3}(i1,i2,i3)];
                                    samp.MQMmtrxC{im1,im2,im3}=[samp.MQMmtrxC{im1,im2,im3},tempC{im1,im2,im3}(i1,i2,i3)];
                                end
                            end
                        end
                    end
                end
            end

            
            clear  i1 i2 i3 im1 im2 im3 Mmtrx_LP tempLP i j tempG tempC

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
    
    eval.u=ones(1,3)*140;eval.sigma2=800;
    tempgg=rand(59,1);
    eval.g=tempgg/sum(tempgg);
    Ou=zeros(1,3);Osigma2=0;
    CRtime=0;
    
    rg=rand(59,1);
    rg=rg/sum(rg);
    eval.g=fmincon(@rand59g,eval.g,[],[],Aeq,Beq,zeros(1,59),ones(1,59),[],optionsSD);
%     (Aeq*eval.g)'
    



    while sum(eval.u-Ou)>0.1||eval.sigma2-Osigma2>1
        [eval.g,fval,exitflag,output]=fmincon(@ComputeL1,eval.g,[],[],Aeq,Beq,zeros(1,59),ones(1,59),[],optionsEP);
%          (Aeq*eval.g)'
        Ou=eval.u;Osigma2=eval.sigma2;  
        part_Em;
        CRtime=CRtime+1;
        if CRtime>15
            break
        end
        fprintf('Simu/Scan %d/%.1f  CRtime=%d--U:(%6.3f  %6.3f  %6.3f)    Sigma2=%6.3f   #emtime=%d\n',simu,scan.point,CRtime,eval.u,eval.sigma2,emtime)
    end
    eval.dr(1)=sum(eval.pai([1,2,3,4]));
    eval.dr(2)=sum(eval.pai([1,2,5,6]));
    eval.dr(3)=sum(eval.g([1,2,5,6,7,12,13,14,19,20,21,26,27,28,33,34,35,40,41,45,46,47,48,55,56]));
    L1=-ComputeL1(eval.g);
    L0=ComputeL0(mean(SM.phen),var(SM.phen));
    LR=-2*(L0-L1);
    
    result=[result,[scan.point;LR;sum(abs(Aeq*eval.g-Beq));0;eval.u';0;eval.sigma2;0;eval.dr']];
    fprintf('**simu/scan %d/%.1f  LR=%.2f   U:(%6.3f  %6.3f  %6.3f)    Sigma2=%6.3f   err=%6.3f\n',simu,scan.point,LR,eval.u,eval.dr,eval.sigma2,sum(abs(Aeq*eval.g-Beq))) 
    
    scan.point=scan.point+SetScan.dens;
    if scan.point>SetScan.end
        break;
    end
    
    clear CRtime L1 L0 LR Osigma2 Ou ou osigma2 etime exitflag fval i j l tempsg tempsum tempu1 tempu2 output freqP tempgg
    
  end



  result1=result(:,result(2,:)==max(result(2,:)));
  fprintf('\n*****simu %d has been done\n',simu);
  fprintf('simu/scan %d/%.1f    Max LR=%.2f  U:%6.3f  %6.3f  %6.3f  Sigma2=%6.3f\n',simu,result1([1,2,5:7,9]));
  fprintf('Dr: %6.3f  %6.3f  %6.3f\n\n',result1(11:13));

  SimuResult(:,simu)=result1;
  save(savename,'SimuResult','samplesize','tval');
  clear result1 result SM Aeqt2 Beqt2
  
end





