


%%%%%%%%%%%%%%%%%%%%%%%%%建立配子指纹

% % gfm{ typical gamete , number , finger print }
% 
% gfm=xlsread('d:\matlab\work\autotetraploidQTL\data\gfm107.xls','sheet1');
% 
% for i=1:length(gfm)
%     gn=gfm(i,2);
%     for j=1:6
%         n=gn-(fix(gn/(10^j))*10^j);
%         fgp{i,1}(7-j)=fix(n/10^(j-1));
%     end
%     fgp{i,2}=[fgp{i,1}(4:6),fgp{i,1}(1:3)];
%     temp11=zeros(5,5);
%     temp12=zeros(5,5);
%     for k1=1:5
%         for k2=(k1+1):6
%             if fgp{i,1}(k1)~=fgp{i,1}(k2)
%                 temp11(k1,(k2-1))=1;
%             end
%             if fgp{i,2}(k1)~=fgp{i,2}(k2)
%                 temp12(k1,(k2-1))=1;
%             end
%         end
%     end
%     fgp{i,3}=temp11;
%     fgp{i,4}=temp12;
% end
% 
% clear gn i j k1 k2 n gn temp11 temp12
% 
% mkindx=[1,1;2,2;3,3;4,4;1,2;1,3;1,4;2,3;2,4;3,4];
% 
% for m1=1:length(mkindx(:,1))
%     for m2=1:length(mkindx(:,1))
%         for m3=1:length(mkindx(:,1))
%             tempg=[mkindx(m1,1),mkindx(m2,1),mkindx(m3,1),mkindx(m1,2),mkindx(m2,2),mkindx(m3,2)];
%             temp2=zeros(5,5);
%             for k1=1:5
%                 for k2=(k1+1):6
%                     if tempg(k1)~=tempg(k2)
%                         temp2(k1,(k2-1))=1;
%                     end
%                 end
%             end
%             for i=1:length(fgp)
%                 if (sum(sum(fgp{i,3}~=temp2))==0)||(sum(sum(fgp{i,4}~=temp2))==0)
% %                     freq(m1,m2,m3)=tg(gfm(i,4))/gfm(i,3);
%                     coef{1}(m1,m2,m3)=gfm(i,4);
%                     coef{2}(m1,m2,m3)=gfm(i,3);
%                 end
%             end
%         end
%     end
% end
% 
% MQMmtrx=coef;
% clear i k1 k2 m1 m2 m3 temp2 tempg coef
% 
% save('fgp');

load fgp

tval.u=[tval.ovau+tval.add(1),tval.ovau+tval.add(2),tval.ovau+tval.add(3),tval.ovau-tval.add(1)-tval.add(2)-tval.add(3),...
        tval.ovau+tval.add(1)+tval.add(2)+tval.dom(1),tval.ovau+tval.add(1)+tval.add(3)+tval.dom(2),tval.ovau-tval.add(2)-tval.add(3)+tval.dom(3)...
        tval.ovau+tval.add(2)+tval.add(3)+tval.dom(4),tval.ovau-tval.add(1)-tval.add(3)+tval.dom(5),tval.ovau-tval.add(1)-tval.add(2)+tval.dom(6)];

% tval.up(1:4)=tval.dr(tval.loc)/4;
% tval.up(5:10)=(1-tval.dr(tval.loc))/6;
% 
% tval.uu=sum(tval.u.*tval.up);
% tval.sigma2=sum(((tval.u-tval.uu).^2.*tval.up))/tval.h2;




%%%设定计算参数的系数

coefx.pai{1}=1:4;    coefx.pai{2}=5:11;    coefx.pai{3}=12:18;
coefx.pai{4}=19:25;  coefx.pai{5}=26:32;   coefx.pai{6}=33:39;
coefx.pai{7}=40:44;  coefx.pai{8}=45:54;   coefx.pai{9}=55:59;

coefx.dr{1}=1:25;                   %A
coefx.dr{2}=[1:11,26:39];           %B
coefx.dr{3}=[1,2,5,6,7,12,13,14,19,20,21,26,27,28,33,34,35,40,41,45,46,47,48,55,56];    %Q

coefx.r{1,1}=[2,4,6,7,10,11,13,14,17,18,20,21,24,25,28,32,33,35,37,41,44,47,48,54,56,59];   %rAQ
coefx.r{1,2}=[3,8,9,15,16,22,23,26,27,34,40,45,46,55];
coefx.r{1,3}=[29,38,42,49,57];
coefx.r{1,4}=[30,31,36,39,43,50,51,52,53,58];

coefx.r{2,1}=[2,4,5,7,9,11,14,18,19,21,23,27,28,31,32,34,35,38,39,41,44,46,48,53,55,57];   %rQB
coefx.r{2,2}=[3,8,10,12,13,20,29,30,36,37,40,45,47,56];
coefx.r{2,3}=[15,24,42,50,59];
coefx.r{2,4}=[16,17,22,25,43,49,51,52,54,58];

coefx.r{3,1}=[5:11,19:25,33:39,52,55:59];           % 1             %rAB
coefx.r{3,2}=[12:18,26:32,40,41,44,49,50,51];       % 1/2
coefx.r{3,3}=[45,46,47,48,53,54];                   % 3/4
coefx.r{3,4}=[42,49,50,43,51,52,49];                % other

Mmtrx{1}= [1,2,2,2,3,3,3,4,4,4;                %Gs
           2,1,2,2,3,4,4,3,3,4;
           2,2,1,2,4,3,4,3,4,3;
           2,2,2,1,4,4,3,4,3,3;
           5,5,6,6,7,8,8,8,8,9;
           5,6,5,6,8,7,8,8,9,8;
           5,6,6,5,8,8,7,9,8,8;
           6,5,5,6,8,8,9,7,8,8;
           6,5,6,5,8,9,8,8,7,8;
           6,6,5,5,9,8,8,8,8,7];


Mmtrx{2}= [ 4,12,12,12,12,12,12,12,12,12;    % x
           12, 4,12,12,12,12,12,12,12,12;
           12,12, 4,12,12,12,12,12,12,12;
           12,12,12, 4,12,12,12,12,12,12;
           12,12,12,12, 6,24,24,24,24, 6;
           12,12,12,12,24, 6,24,24, 6,24;
           12,12,12,12,24,24, 6, 6,24,24;
           12,12,12,12,24,24, 6, 6,24,24;
           12,12,12,12,24, 6,24,24, 6,24;
           12,12,12,12, 6,24,24,24,24, 6];
%%%染色体

scan.chr=0;tempsum=0;
for i=1:length(tval.map)
    tempsum=tempsum+tval.map(i);
    scan.chr=[scan.chr,tempsum];
end
scan.chr=[scan.chr,0];
clear tempsum i
 
       
       
