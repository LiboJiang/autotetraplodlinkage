


library(mappoly)

tetra.solcap <- read_geno_csv(file.in  = "../Data/tetra_1_solcap.csv", ploidy = 4)

all.mrk <- make_seq_mappoly(tetra.solcap, 'all')
red.mrk <- elim_redundant(all.mrk)
unique.mrks <- make_seq_mappoly(red.mrk)
counts <- cache_counts_twopt(unique.mrks, cached = TRUE)

all.pairs <- est_pairwise_rf(input.seq = unique.mrks,
                             count.cache = counts,
                             n.clusters = 1, 
                             verbose=TRUE)


seq1_map <-est_rf_hmm_sequential(input.seq = unique.mrks,
                           start.set = 10,
                           thres.twopt = 10, 
                           thres.hmm = 10,
                           extend.tail = 200,
                           info.tail = TRUE, 
                           twopt = all.pairs,
                           sub.map.size.diff.limit = 8, 
                           phase.number.limit = 10,
                           reestimate.single.ph.configuration = TRUE,
                           tol = 10e-3,
                           tol.final = 10e-4)


seq1_map1 <- est_full_hmm_with_global_error(input.map = seq1_map,error = 0.05, 
                                  tol = 10e-4, verbose = TRUE)


plot_map_list(seq1_map1)

d <- extract_map(seq1_map1)
z <- data.frame(mrk = get(seq1_map1$info$data.name)$mrk.names[seq1_map1$maps[[1]]$seq.num], 
                         LG = rep("lg1",length(d)), pos = d)


colnames(z) <- c("snp_name","LG","pos")
write.csv(z,file="HMM_LG1.csv")


