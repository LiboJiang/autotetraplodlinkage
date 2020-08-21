
%%% figer our Linkage Phase & Pai
%%% for phase 23,34,45,56

function RetureValue = ConfirmLP(PHcase)

global Mmtrx samp Lphase samplesize SM

for mk1=1:Lphase(PHcase).dms(1)
    for mk2=1:Lphase(PHcase).dms(2)
        MtrxPaiP{mk1,mk2}=[];
        MtrxCoef{mk1,mk2}=[];
        for i1=1:length(Lphase(PHcase).mtc{1,mk1})
            for i2=1:length(Lphase(PHcase).mtc{2,mk2})
                MtrxPaiP{mk1,mk2}=[MtrxPaiP{mk1,mk2},Mmtrx{1}(Lphase(PHcase).mtc{1,mk1}(i1),Lphase(PHcase).mtc{2,mk2}(i2))];
                MtrxCoef{mk1,mk2}=[MtrxCoef{mk1,mk2},Mmtrx{2}(Lphase(PHcase).mtc{1,mk1}(i1),Lphase(PHcase).mtc{2,mk2}(i2))];
            end
        end
    end
end



opai=zeros(9,1);
epai=zeros(9,1)+1/9;
egtime=0;

while sum(abs(opai-epai)>0.000001)>1
opai=epai;

for ipai=1:length(epai)
    tempss=0;
    for mk1=1:Lphase(PHcase).dms(1)
        for mk2=1:Lphase(PHcase).dms(2)
            if (sum(MtrxPaiP{mk1,mk2}==ipai)>0)&&(samp.Nmtrx(mk1,mk2)>0)
                tempss=tempss+sum(   epai(ipai)./MtrxCoef{mk1,mk2}(MtrxPaiP{mk1,mk2}==ipai)   )     /    sum(   epai(MtrxPaiP{mk1,mk2})'./MtrxCoef{mk1,mk2}   )*samp.Nmtrx(mk1,mk2);   
            end
        end
    end
    epai(ipai)=tempss/samplesize;    
end
egtime=egtime+1;
end


LphasePar(PHcase).epai=epai';

output.epai=LphasePar(PHcase).epai;
output.LinkagePhase=Lphase(PHcase).case;
output.mtc=Lphase(PHcase).mtc;
output.Mpai{1}=MtrxPaiP;
output.Mpai{2}=MtrxCoef;


fprintf('LinkagePhase:[%s]  epai: %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f %4.3f \n',Lphase(PHcase).case,epai(1),epai(2),epai(3),epai(4),epai(5),epai(6),epai(7),epai(8),epai(9));

clear efreq epai



RetureValue=output;








