       



for l=1:length(samp.PHmtrx{i,j})
            xx0=log(normpdf(samp.PHmtrx{i,j}(1),tempu,sqrt(tempsigma2)));
end


for l=1:length(samp.PHmtrx{i,j})
            xx1=log(sum(    eval.g(MQMmtrx{1}(i,:,j))'./MQMmtrx{2}(i,:,j)/(eval.pai(Mmtrx{1}(i,j))/Mmtrx{2}(i,j)).* normpdf(samp.PHmtrx{i,j}(1),eval.u,sqrt(eval.sigma2))    ));                                  
end