function value = ComputeL0(tempu,tempsigma2)
% for cpt L0
global samp mdm

tempL0=zeros(10,10);

for i=1:mdm(1)
    for j=1:mdm(2)
        for l=1:length(samp.PHmtrx{i,j})
            tempL0(i,j)=tempL0(i,j)+log(normpdf(samp.PHmtrx{i,j}(l),tempu,sqrt(tempsigma2)));
        end
    end
end

logL0=sum(sum(tempL0));
value=logL0;

end

