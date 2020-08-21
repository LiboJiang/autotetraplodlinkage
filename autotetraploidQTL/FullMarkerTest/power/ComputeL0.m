function value = ComputeL0(tempu,tempsigma2)
% for cpt L0
global samp

tempL0=zeros(10,10);

for i=1:10
    for j=1:10
        for l=1:length(samp.PHmtrx{i,j})
            tempL0(i,j)=tempL0(i,j)+log(normpdf(samp.PHmtrx{i,j}(l),tempu,sqrt(tempsigma2)));
        end
    end
end

logL0=sum(sum(tempL0));
value=logL0;

end

