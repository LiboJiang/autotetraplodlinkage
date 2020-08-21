
%%%%%for scan QTl


clc;
clear;
clear global;

warning off all

global Mmtrx samp Lphase eval rg
% loaddataSQ;
load SQdata
optionsEP=optimset('Algorithm','interior-point','DiffMinChange',1e-12,'MaxIter',500,'MaxFunEvals',5000,'Display','off');
optionsSD=optimset('Algorithm','interior-point', 'Display','off');%,'MaxFunEvals',5000);
rand('state',sum(100*clock));

%%%% Set parameter
SetScan.start=0;
SetScan.end=10000;
SetScan.dens=0.5;

% 1000
qphase{1}=1;
qphase{2}=[5,6,7];
qphase{3}=[2,3,4,8,9,10];

scan.point=SetScan.start;
scan.mk=0;scan.sideL=0;scan.sideR=-1;
result=[];

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
        samp.size=0;samp.phenotype=[];samp.marker=[];samp.mtrx=zeros(2,2);
        samp.PhenMtrx{1,1}=[];samp.PhenMtrx{1,2}=[];samp.PhenMtrx{2,1}=[];samp.PhenMtrx{2,2}=[];
        for i=1:(length(ogdata.marker(:,1))-1)
             if (ogdata.marker(i+1,scan.mk)>3)||(ogdata.marker(i+1,scan.mk+1)>3)||(isnan(ogdata.phenotype(i+1,2)));
                 continue;
             end
             samp.size=samp.size+1;
             samp.phenotype=[samp.phenotype;ogdata.phenotype(i+1,2)];
             samp.marker(samp.size,:)=(ogdata.marker(i+1,scan.mk:scan.mk+1));
             samp.mtrx(-samp.marker(samp.size,1)+2,-samp.marker(samp.size,2)+2)=samp.mtrx(-samp.marker(samp.size,1)+2,-samp.marker(samp.size,2)+2)+1;
             samp.PhenMtrx{-samp.marker(samp.size,1)+2,-samp.marker(samp.size,2)+2}=[samp.PhenMtrx{-samp.marker(samp.size,1)+2,-samp.marker(samp.size,2)+2},ogdata.phenotype(i+1,2)];
        end
        
        %确定连锁相,计算pai
        tempLP=ConfirmLP;
        eval.pai=tempLP.epai;
        eval.LinkagePhase=tempLP.LinkagePhase;
        Mmtrx_LP=zeros(2,2);
        for i1=1:2
            for i2=1:2
                Mmtrx_LP(i1,i2)=sum(eval.pai(tempLP.Mpai{1}{i1,i2})./tempLP.Mpai{2}{i1,i2}');
            end
        end
        samp.Pai_Mmtrx=Mmtrx_LP;
        

        %压缩samp.MQM矩阵      
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
    Beq=[eval.pai;eval.r';1];
    
   
%     test=20;
%     testm=[];
%     
%     
%     for i=1:9
%         if Beq(i)<0.00001
%         tempg(coefx.pai{i})=0;
%         end
%     end
%     
% %     for i=1:2
% %         if Beq(9+i)<0.00001
% %         tempg(coefx.r{i,1})=0;tempg(coefx.r{i,2})=0;tempg(coefx.r{i,3})=0;tempg(coefx.r{i,4})=0;
% %         end
% %     end
% %     

   
    rg=rand(59,1);
    rg=rg/sum(rg);
    tempg=rand(59,1);
    tempg=fmincon(@rand59g,tempg,[],[],Aeq,Beq,zeros(1,59),ones(1,59),[],optionsSD);
%     Aeq*tempg-Beq

   

    %%% est para
    eval.g=tempg;
    eval.u=ones(1,3)*150;eval.sigma2=700;eval.g=tempg;
    Ou=zeros(1,3);Osigma2=0;
    CRtime=0;
    

    while sum(eval.u-Ou)>0.0001||eval.sigma2-Osigma2>0.1
        [eval.g,fval,exitflag,output]=fmincon(@ComputeL1,eval.g,[],[],Aeq,Beq,zeros(1,59),ones(1,59),[],optionsEP);
%         (Aeq*tempg-Beq)'
        Ou=eval.u;Osigma2=eval.sigma2;
        part_EM;
        CRtime=CRtime+1;
        fprintf('Scan %.1f  CRtime=%d--U:(%6.3f  %6.3f  %6.3f)    Sigma2=%6.3f   #emtime=%d\n',scan.point,CRtime,eval.u,eval.sigma2,emtime)
        if CRtime>15
            break
        end
    end
    
    eval.dr(1)=sum(eval.pai([1,2,3,4]));
    eval.dr(2)=sum(eval.pai([1,2,5,6]));
    eval.dr(3)=sum(eval.g([1,2,5,6,7,12,13,14,19,20,21,26,27,28,33,34,35,40,41,45,46,47,48,55,56]));
    
    L1=-ComputeL1(eval.g);
    L0=ComputeL0(mean(samp.phenotype),var(samp.phenotype));
    LR=-2*(L0-L1);
    result=[result,[scan.point;LR;sum(abs(Aeq*eval.g-Beq));0;eval.u';0;eval.sigma2;0;eval.dr';eval.LinkagePhase]];

    fprintf('**Scan %.1f  LR=%.2f   U:(%6.3f  %6.3f  %6.3f)    Sigma2=%6.3f  Phase=%d  err=%6.3f\n\n',scan.point,LR,eval.u,eval.sigma2,eval.LinkagePhase,sum(abs(Aeq*eval.g-Beq))) 
    

    scan.point=scan.point+SetScan.dens;
    if scan.point>SetScan.end
        break;
    end
    
    clear CRtime L1 L0 LR Osigma2 Ou ou osigma2 etime exitflag fval i j l tempsg tempsum tempu1 tempu2 output freqP tempg emtime
    
    save('ScanResult','result');
    
end
  

  result1=result(:,result(2,:)==max(result(2,:)));
  fprintf('\nScan has been done\n');
  fprintf('Scan %.1f    Max LR=%.2f  U:%6.3f  %6.3f  %6.3f    Sigma2=%6.3f     Phase=%d\n',result1([1,2,5:7,9,14]));
  fprintf('Dr: %6.3f  %6.3f  %6.3f\n',result1(11:13));


  

    
    
    
    
    
    
