function value = ComputeL1(og)
% for cpt L1

global eval Mmtrx samp MQMmtrx
eval.g=og;
tempL1=zeros(10,10);

for i=1:10
    for j=1:10
        for l=1:length(samp.PHmtrx{i,j})
            tempL1(i,j)=tempL1(i,j)+log(sum(    eval.g(MQMmtrx{1}(i,:,j))'./MQMmtrx{2}(i,:,j)/(eval.pai(Mmtrx{1}(i,j))/Mmtrx{2}(i,j)).* normpdf(samp.PHmtrx{i,j}(l),eval.u,sqrt(eval.sigma2))));
        end
    end
end

logL1=sum(sum(tempL1));
value=-logL1;

end

