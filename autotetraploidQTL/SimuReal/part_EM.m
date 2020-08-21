


%%%% for em cpt u & sigma2



ou=zeros(1,3);osigma2=0;
emtime=0;


while sum(eval.u-ou)>0.01||eval.sigma2-osigma2>1

tempu1=zeros(1,3);
tempu2=zeros(1,3);
tempsg=0;

for i=1:2
    for j=1:2
        for l=1:length(samp.PhenMtrx{i,j})
            tempsum=sum(  [sum(eval.g(samp.MQMmtrxG{i,1,j})'./samp.MQMmtrxC{i,1,j}),...
                           sum(eval.g(samp.MQMmtrxG{i,2,j})'./samp.MQMmtrxC{i,2,j}),...
                           sum(eval.g(samp.MQMmtrxG{i,3,j})'./samp.MQMmtrxC{i,3,j})]  .*normpdf(samp.PhenMtrx{i,j}(l),eval.u,sqrt(eval.sigma2)));
            freqP{i,j}(l,:)=[sum(eval.g(samp.MQMmtrxG{i,1,j})'./samp.MQMmtrxC{i,1,j}),...
                             sum(eval.g(samp.MQMmtrxG{i,2,j})'./samp.MQMmtrxC{i,2,j}),...
                             sum(eval.g(samp.MQMmtrxG{i,3,j})'./samp.MQMmtrxC{i,3,j})] .*normpdf(samp.PhenMtrx{i,j}(l),eval.u,sqrt(eval.sigma2))/tempsum;
         
            tempu1=tempu1+freqP{i,j}(l,:).*samp.PhenMtrx{i,j}(l);
            tempu2=tempu2+freqP{i,j}(l,:);
            tempsg=tempsg+sum(freqP{i,j}(l,:).*(samp.PhenMtrx{i,j}(l)-eval.u).^2);
        end
    end
end

ou=eval.u;
osigma2=eval.sigma2;

eval.u=tempu1./tempu2;
eval.sigma2=tempsg/samplesize;

emtime=emtime+1;

end



