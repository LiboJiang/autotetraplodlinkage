setwd("E:/Pro1/AutoL1/HMM_s1")

library(mappoly)


#5

all.mrk5 <- make_seq_mappoly(tetra.solcap, 'seq5')
red.mrk5 <- elim_redundant(all.mrk5)
unique.mrks5 <- make_seq_mappoly(red.mrk5)
counts5 <- cache_counts_twopt(unique.mrks5, cached = TRUE)

all.pairs5 <- est_pairwise_rf(input.seq = unique.mrks5,
                             count.cache = counts5,
                             verbose=TRUE)


seq1_map5 <-est_rf_hmm_sequential(input.seq = unique.mrks5,
                                 start.set = 10,
                                 thres.twopt = 10, 
                                 thres.hmm = 10,
                                 extend.tail = 200,
                                 info.tail = TRUE, 
                                 twopt = all.pairs5,
                                 sub.map.size.diff.limit = 8, 
                                 phase.number.limit = 20,
                                 reestimate.single.ph.configuration = F,
                                 tol = 10e-3,
                                 tol.final = 10e-4)

plot_map_list(seq1_map5)

d5 <- extract_map(seq1_map5)
z5 <- data.frame(mrk = get(seq1_map5$info$data.name)$mrk.names[seq1_map5$maps[[1]]$seq.num], 
                LG = rep("lg5",length(d5)), pos = d5)


write.csv(z5,file="LG5_Mappoly.csv")


ph5_p <- ph_list_to_matrix(seq1_map5$maps[[1]]$seq.ph$P,4)
ph5_q <- ph_list_to_matrix(seq1_map5$maps[[1]]$seq.ph$Q,4)

ph5 <- cbind(ph5_p,ph5_q)
write.csv(ph5,file="phase5.csv")


#9

all.mrk9 <- make_seq_mappoly(tetra.solcap, 'seq9')
red.mrk9 <- elim_redundant(all.mrk9)
unique.mrks9 <- make_seq_mappoly(red.mrk9)
counts9 <- cache_counts_twopt(unique.mrks9, cached = TRUE)

all.pairs9 <- est_pairwise_rf(input.seq = unique.mrks9,
                              count.cache = counts9,
                              verbose=TRUE)


seq1_map9 <-est_rf_hmm_sequential(input.seq = unique.mrks9,
                                  start.set = 10,
                                  thres.twopt = 10, 
                                  thres.hmm = 10,
                                  extend.tail = 200,
                                  info.tail = TRUE, 
                                  twopt = all.pairs9,
                                  sub.map.size.diff.limit = 8, 
                                  phase.number.limit = 10,
                                  reestimate.single.ph.configuration = F,
                                  tol = 10e-3,
                                  tol.final = 10e-4)

plot_map_list(seq1_map9)

d9 <- extract_map(seq1_map9)
z9 <- data.frame(mrk = get(seq1_map9$info$data.name)$mrk.names[seq1_map9$maps[[1]]$seq.num], 
                 LG = rep("lg9",length(d9)), pos = d9)



write.csv(z9,file="LG9_Mappoly.csv")

ph9_p <- ph_list_to_matrix(seq1_map9$maps[[1]]$seq.ph$P,4)
ph9_q <- ph_list_to_matrix(seq1_map9$maps[[1]]$seq.ph$Q,4)

ph9 <- cbind(ph9_p,ph9_q)
write.csv(ph9,file="phase9.csv")
