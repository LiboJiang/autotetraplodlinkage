function value = ComputeL0(tempu,tempsigma2)
% for cpt L0
global samp

tempL0=zeros(10,10);

for i=1:2
    for j=1:2
        if isempty(samp.PhenMtrx{i,j})
            continue;
        end
        for l=1:length(samp.PhenMtrx{i,j})
            tempL0(i,j)=tempL0(i,j)+log(normpdf(samp.PhenMtrx{i,j}(l),tempu,sqrt(tempsigma2)));
        end
    end
end

logL0=sum(sum(tempL0));
value=logL0;

end

