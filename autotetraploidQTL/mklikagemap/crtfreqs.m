

%%%
%%% 产生不同情况样本矩阵
%%%

dms=[2,2,2];


%%% case 11 
%%% 1 0 0 0
%%% 1 0 0 0
%%% 1 0 0 0

mtc{1,1}=[1,5,6,7];         mtc{1,2}=[2,3,4,8,9,10];  
mtc{2,1}=[1,5,6,7];         mtc{2,2}=[2,3,4,8,9,10];     
mtc{3,1}=[1,5,6,7];         mtc{3,2}=[2,3,4,8,9,10];       

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs11=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 12 
%%% 1 0 0 0
%%% 1 0 0 0
%%% 0 1 0 0

mtc{1,1}=[1,5,6,7];         mtc{1,2}=[2,3,4,8,9,10];  
mtc{2,1}=[1,5,6,7];         mtc{2,2}=[2,3,4,8,9,10];     
mtc{3,1}=[2,5,8,9];         mtc{3,2}=[1,3,4,6,7,10];       

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs12=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 13 
%%% 1 0 0 0
%%% 0 1 0 0
%%% 1 0 0 0

mtc{1,1}=[1,5,6,7];         mtc{1,2}=[2,3,4,8,9,10];  
mtc{2,1}=[2,5,8,9];         mtc{2,2}=[1,3,4,6,7,10];     
mtc{3,1}=[1,5,6,7];         mtc{3,2}=[2,3,4,8,9,10];       

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs13=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 14 
%%% 1 0 0 0
%%% 0 1 0 0
%%% 0 1 0 0

mtc{1,1}=[1,5,6,7];         mtc{1,2}=[2,3,4,8,9,10];  
mtc{2,1}=[2,5,8,9];         mtc{2,2}=[1,3,4,6,7,10];     
mtc{3,1}=[2,5,8,9];         mtc{3,2}=[1,3,4,6,7,10];       

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs14=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 21 
%%% 1 1 0 0
%%% 1 1 0 0
%%% 1 1 0 0

mtc{1,1}=[1,2,5,6,7,8,9];         mtc{1,2}=[3,4,10];  
mtc{2,1}=[1,2,5,6,7,8,9];         mtc{2,2}=[3,4,10];    
mtc{3,1}=[1,2,5,6,7,8,9];         mtc{3,2}=[3,4,10];       

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs21=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 22 
%%% 1 1 0 0
%%% 1 1 0 0
%%% 1 0 1 0

mtc{1,1}=[1,2,5,6,7,8,9];         mtc{1,2}=[3,4,10];  
mtc{2,1}=[1,2,5,6,7,8,9];         mtc{2,2}=[3,4,10];    
mtc{3,1}=[1,3,5,6,7,8,10];        mtc{3,2}=[2,4,9];          

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs22=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 23 
%%% 1 1 0 0
%%% 1 0 1 0
%%% 1 1 0 0


mtc{1,1}=[1,2,5,6,7,8,9];         mtc{1,2}=[3,4,10];  
mtc{2,1}=[1,3,5,6,7,8,10];        mtc{2,2}=[2,4,9];    
mtc{3,1}=[1,2,5,6,7,8,9];         mtc{3,2}=[3,4,10];       

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs23=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 24
%%% 1 1 0 0
%%% 1 0 1 0
%%% 1 0 1 0

mtc{1,1}=[1,2,5,6,7,8,9];         mtc{1,2}=[3,4,10];  
mtc{2,1}=[1,3,5,6,7,8,10];        mtc{2,2}=[2,4,9];    
mtc{3,1}=[1,3,5,6,7,8,10];        mtc{3,2}=[2,4,9];      

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs24=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 25 
%%% 1 1 0 0
%%% 1 0 1 0
%%% 1 0 0 1

mtc{1,1}=[1,2,5,6,7,8,9];         mtc{1,2}=[3,4,10];  
mtc{2,1}=[1,3,5,6,7,8,10];        mtc{2,2}=[2,4,9];    
mtc{3,1}=[1,4,5,6,7,9,10];        mtc{3,2}=[2,3,8];    

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs25=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 26 
%%% 1 1 0 0
%%% 1 1 0 0
%%% 0 0 1 1

mtc{1,1}=[1,2,5,6,7,8,9];         mtc{1,2}=[3,4,10]; 
mtc{2,1}=[1,2,5,6,7,8,9];         mtc{2,2}=[3,4,10];    
mtc{3,1}=[3,4,6,7,8,9,10];        mtc{3,2}=[1,2,5];       

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs26=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 27 
%%% 1 1 0 0
%%% 1 0 1 1
%%% 1 1 0 0

mtc{1,1}=[1,2,5,6,7,8,9];         mtc{1,2}=[3,4,10];   
mtc{2,1}=[1,3,4,5,6,7,8,9,10];    mtc{2,2}=2;     
mtc{3,1}=[1,2,5,6,7,8,9];         mtc{3,2}=[3,4,10];        

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs27=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 28 
%%% 1 1 0 0
%%% 0 0 1 1
%%% 0 0 1 1

mtc{1,1}=[1,2,5,6,7,8,9];         mtc{1,2}=[3,4,10]; 
mtc{2,1}=[3,4,6,7,8,9,10];        mtc{2,2}=[1,2,5];   
mtc{3,1}=[3,4,6,7,8,9,10];        mtc{3,2}=[1,2,5];       
      

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs28=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 31 
%%% 1 1 1 0
%%% 1 1 1 0
%%% 1 1 1 0

mtc{1,1}=[1,2,3,5,6,7,8,9,10];         mtc{1,2}=4;  
mtc{2,1}=[1,2,3,5,6,7,8,9,10];         mtc{2,2}=4;     
mtc{3,1}=[1,2,3,5,6,7,8,9,10];         mtc{3,2}=4;       

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs31=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 32 
%%% 1 1 1 0
%%% 1 1 1 0
%%% 1 1 0 1

mtc{1,1}=[1,2,3,5,6,7,8,9,10];         mtc{1,2}=4;  
mtc{2,1}=[1,2,3,5,6,7,8,9,10];         mtc{2,2}=4;     
mtc{3,1}=[1,2,4,5,6,7,8,9,10];         mtc{3,2}=3;      

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs32=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 33 
%%% 1 1 1 0
%%% 1 1 0 1
%%% 1 1 1 0

mtc{1,1}=[1,2,3,5,6,7,8,9,10];         mtc{1,2}=4;  
mtc{2,1}=[1,2,4,5,6,7,8,9,10];         mtc{2,2}=3;     
mtc{3,1}=[1,2,3,5,6,7,8,9,10];         mtc{3,2}=4;        

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs33=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3

%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%

%%% case 34 
%%% 1 1 1 0
%%% 1 1 0 1
%%% 1 1 0 1

mtc{1,1}=[1,2,3,5,6,7,8,9,10];         mtc{1,2}=4;  
mtc{2,1}=[1,2,4,5,6,7,8,9,10];         mtc{2,2}=3;     
mtc{3,1}=[1,2,4,5,6,7,8,9,10];         mtc{3,2}=3;      

for mk1=1:dms(1)
    for mk2=1:dms(2)
        for mk3=1:dms(3)
            coefs{1,mk1,mk2,mk3}=[];
            coefs{2,mk1,mk2,mk3}=[];
            for i1=1:length(mtc{1,mk1})
                for i2=1:length(mtc{2,mk2})
                    for i3=1:length(mtc{3,mk3})
                        coefs{1,mk1,mk2,mk3}=[coefs{1,mk1,mk2,mk3},coef{1}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                        coefs{2,mk1,mk2,mk3}=[coefs{2,mk1,mk2,mk3},coef{2}(mtc{1,mk1}(i1),mtc{2,mk2}(i2),mtc{3,mk3}(i3))];
                    end
                end
            end
        end
    end
end
coefs34=coefs;
clear mtc mk1 mk2 mk3 i1 i2 i3






















