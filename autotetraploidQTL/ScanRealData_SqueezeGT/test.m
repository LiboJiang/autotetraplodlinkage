
sum(sum(eval.pai(Mmtrx{1})./Mmtrx{2}))


ss=0;

for iii=1:2
    for jjj=1:2
        ss=ss+sum(sum(sum(eval.g(samp.MQMmtrxG{iii,jjj})./samp.MQMmtrxC{iii,jjj})));
    end
end
        