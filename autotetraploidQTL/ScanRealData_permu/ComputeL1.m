function value = ComputeL1(og)
% for cpt L1

global eval samp 

tempL1=zeros(2,2);
og=og';

for i=1:2
    for j=1:2
        if isempty(samp.PhenMtrx{i,j})
            continue;
        end
        for l=1:length(samp.PhenMtrx{i,j})
            tempL1(i,j)=tempL1(i,j)+log(...
            sum(  [sum(og(samp.MQMmtrxG{i,1,j})./samp.MQMmtrxC{i,1,j}),...
                   sum(og(samp.MQMmtrxG{i,2,j})./samp.MQMmtrxC{i,2,j}),...
                   sum(og(samp.MQMmtrxG{i,3,j})./samp.MQMmtrxC{i,3,j})]...
            /samp.Pai_Mmtrx(i,j)...
            .* normpdf(samp.PhenMtrx{i,j}(l),eval.u,sqrt(eval.sigma2))...
            )...
            );
        end
    end
end

logL1=sum(sum(tempL1));
value=-logL1;

end

