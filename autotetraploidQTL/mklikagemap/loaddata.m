
%%%%%%%%%%%%new real data%%%%%%%%%%%%%


clear;
clc;


[ogdatamk,ogmkn]=xlsread('d:\matlab\work\autotetraploid\data\jointmapdata.xls','k5test'); 
mkn=ogmkn(2:length(ogmkn(:,1)));
datamk=ogdatamk';


% % sn=0;
% % for i=1:length(ogdatamk(1,:))
% %     if ogdatamk(1,i)==0&&ogdatamk(2,i)==1&&sum(ogdatamk(:,i)==9)<20;
% %         sn=sn+1;
% %         datamk(:,sn)=ogdatamk(:,i);mkname(sn)=ogmkname(i);
% %     end
% % end

clear ogdatamk ogmkn
        



%%%%%%%%%%%%%%%%%%%%%%%%%%建立配子指纹

% gfm{ typical gamete , number , finger print }

gfm=xlsread('d:\matlab\work\autotetraploid\data\gfm107.xls','sheet1');

for i=1:length(gfm)
    gn=gfm(i,2);
    for j=1:6
        n=gn-(fix(gn/(10^j))*10^j);
        fgp{i,1}(7-j)=fix(n/10^(j-1));
    end
    fgp{i,2}=[fgp{i,1}(4:6),fgp{i,1}(1:3)];
    temp11=zeros(5,5);
    temp12=zeros(5,5);
    for k1=1:5
        for k2=(k1+1):6
            if fgp{i,1}(k1)~=fgp{i,1}(k2)
                temp11(k1,(k2-1))=1;
            end
            if fgp{i,2}(k1)~=fgp{i,2}(k2)
                temp12(k1,(k2-1))=1;
            end
        end
    end
    fgp{i,3}=temp11;
    fgp{i,4}=temp12;
end

clear gn i j k1 k2 n gn temp11 temp12


%%%设定计算参数的系数


ax=1:25;
bx=[1:11,26:39];
cx=[1,2,5,6,7,12,13,14,19,20,21,26,27,28,33,34,35,40,41,45,46,47,48,55,56];

r1x{1}=[2,4,6,7,10,11,13,14,17,18,20,21,24,25,28,32,33,35,37,41,44,47,48,54,56,59];
r1x{2}=[3,8,9,15,16,22,23,26,27,34,40,45,46,55];
r1x{3}=[29,38,42,49,57];
r1x{4}=[30,31,36,39,43,50,51,52,53,58];

r2x{1}=[2,4,5,7,9,11,14,18,19,21,23,27,28,31,32,34,35,38,39,41,44,46,48,53,55,57];
r2x{2}=[3,8,10,12,13,20,29,30,36,37,40,45,47,56];
r2x{3}=[15,24,42,50,59];
r2x{4}=[16,17,22,25,43,49,51,52,54,58];

r3x{1}=[5:11,19:25,33:39,52,55:59];           % 1
r3x{2}=[12:18,26:32,40,41,44,49,50,51];       % 1/2
r3x{3}=[45,46,47,48,53,54];                   % 3/4
r3x{4}=[42,49,50,43,51,52,49];                % other



%%%%%%%%%%%%%%%%%%%%%%生成频率矩阵
mkindx=[1,1;2,2;3,3;4,4;1,2;1,3;1,4;2,3;2,4;3,4];

for m1=1:length(mkindx(:,1))
    for m2=1:length(mkindx(:,1))
        for m3=1:length(mkindx(:,1))
            tempg=[mkindx(m1,1),mkindx(m2,1),mkindx(m3,1),mkindx(m1,2),mkindx(m2,2),mkindx(m3,2)];
            temp2=zeros(5,5);
            for k1=1:5
                for k2=(k1+1):6
                    if tempg(k1)~=tempg(k2)
                        temp2(k1,(k2-1))=1;
                    end
                end
            end
            for i=1:length(fgp)
                if (sum(sum(fgp{i,3}~=temp2))==0)||(sum(sum(fgp{i,4}~=temp2))==0)
                    coef{1}(m1,m2,m3)=gfm(i,4);
                    coef{2}(m1,m2,m3)=gfm(i,3);
                end
            end
        end
    end
end
                    
    

clear i k1 k2 m1 m2 m3 mkindx temp11 temp12 tempg temp2 test gfm fgp 

save('k5Tdatamk');



