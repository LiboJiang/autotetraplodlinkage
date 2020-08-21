
function OutMarker = SimuMarker2p(Pmarker,dr,d)

 global rg Mmtrx

 optionsSD=optimset('Algorithm','interior-point', 'Display','off');%,'MaxFunEvals',5000);


 tempr=3/4*(1-exp(-4*( d /100)/3));
 tempphi=tempr^2/(9*(1-tempr)^2+tempr^2);
 temppsi=tempr/(3-2*tempr);

 Aeqt2(1,:)=[1,1,1,1,0,0,0,0,0];
 Aeqt2(2,:)=[1,1,0,0,1,1,0,0,0];
 Aeqt2(3,:)=[0,1,1/2,1,1/2,1,tempphi,(1+temppsi)/2,1];
 Aeqt2(4,:)=1;
 Beqt2=[dr(1);dr(2);tempr;1];

 rg=rand(1,9);
 rg=rg/sum(rg);
 tempg=rand(1,9);
 tempg=fmincon(@rand9g,tempg,[],[],Aeqt2,Beqt2,zeros(1,9),ones(1,9),[],optionsSD);
 tempfreq=tempg(Mmtrx{1})./Mmtrx{2};


 for i=1:length(Pmarker)
    tempj=rand(1)*sum(tempfreq(Pmarker(i),:));tempsum=0;
     for j=1:10
        tempsum=tempsum+tempfreq(Pmarker(i),j);
         if tempj<tempsum
            OutMarker(i)=j;tempj=9999;
         end
     end
 end
            
end




















