% %%%
% %%% 计算参数
% %%%


%%%%%%%%%%%%%%%%%%%em估算Gs

dms=size(samp);
og=zeros(59,1);
eg=zeros(59,1)+1/59;
% egtime=0;

tt=zeros(59,1);

while sum(abs(og-eg)>0.000001)>1
og=eg;

for ig=1:59
    tempss=0;
    for mk1=1:dms(1)
        for mk2=1:dms(2)
            for mk3=1:dms(3)
                if (sum(coefs{1,mk1,mk2,mk3}==ig)>0)&&(samp(mk1,mk2,mk3)>0)
                    tempss=tempss+sum(eg(ig)./coefs{2,mk1,mk2,mk3}(coefs{1,mk1,mk2,mk3}==ig))/sum(eg(coefs{1,mk1,mk2,mk3})'./coefs{2,mk1,mk2,mk3})*samp(mk1,mk2,mk3);   
                end
            end
        end
    end
    eg(ig)=tempss/samplesize;    
end

% egtime=egtime+1

end




clear tempss mk1 mk2 mk3 ig og egtime
                

%%%%%%%%%%%%%%%%%%em迭代估算参数


ophi1=0;ophi2=0;opsi1=0;opsi2=0;         ephi1=0.2;ephi2=0.2;epsi1=0.2;epsi2=0.2;
or(1)=0;or(2)=0;                         er(1)=0.2;er(2)=0.2;

eptime=0;

while abs(ophi1-ephi1)>0.000001||abs(ophi2-ephi2)>0.000001||abs(opsi1-epsi1)>0.000001||abs(opsi2-epsi2)>0.000001||abs(or(1)-er(1))>0.000001||abs(or(2)-er(2))>0.000001
    or=er;
    er(1)=sum(eg(r1x{1}))+sum(eg(r1x{2}))/2+sum(eg(r1x{3}))*ephi1+sum(eg(r1x{4}))*(1+epsi1)/2;
    er(2)=sum(eg(r2x{1}))+sum(eg(r2x{2}))/2+sum(eg(r2x{3}))*ephi2+sum(eg(r2x{4}))*(1+epsi2)/2;
    
    ophi1=ephi1;ophi2=ephi2;opsi1=epsi1;opsi2=epsi2;
    ephi1=er(1)^2/(9-18*er(1)+10*er(1)^2);
    ephi2=er(2)^2/(9-18*er(2)+10*er(2)^2);
    epsi1=er(1)/(3-2*er(1));
    epsi2=er(2)/(3-2*er(2));

   
end

ea(1)=sum(eg(ax));ea(2)=sum(eg(bx));ea(3)=sum(eg(cx));
er(3)=sum(eg(r3x{1}))+sum(eg(r3x{2}))/2+sum(eg(r3x{3}))*3/4+(eg(42)*2+eg(49))*ephi1/2+(eg(42)*2+eg(50))*ephi2/2-eg(42)*ephi1*ephi2*2+...
        (eg(43)*2+eg(50)+eg(51)-eg(52))*epsi1/2-eg(50)*ephi2*epsi1+(eg(43)*2+eg(49)+eg(51)-eg(52))*epsi2/2-eg(49)*ephi1*epsi2-(2*eg(43)+eg(51)-eg(52))*epsi2*epsi1; 
   


for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            efreqs(mk1,mk2,mk3)=sum(eg(coefs{1,mk1,mk2,mk3})'./coefs{2,mk1,mk2,mk3});
        end
    end
end
lh=sum(sum(sum(samp.*log(efreqs))));


clear eptime ophi1 ophi2 opsi1 opsi2 or mk1 mk2 mk3











































































