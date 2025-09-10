
log(p)+log(q)

sqrt(800)


n=800
p=70
q=120
nu=0.3
lambda =80
path = paste0("datarepo/n=", n, "p=", p, "q=", q, "nu=", nu)
post_bias <- as.matrix(read.csv(file.path(path, "postselection_bias.csv"),
                                header = TRUE,
                                na.strings = c("NA", "NaN")))

post_len <- as.matrix(read.csv(file.path(path, "postselection_len.csv"),
                               header = TRUE,
                               na.strings = c("NA", "NaN")))

post_count <- as.matrix(read.csv(file.path(path, "postselection_count.csv"),
                                 header = TRUE,
                                 na.strings = c("NA", "NaN")))

post_rej <- as.matrix(read.csv(file.path(path, "postselection_rej.csv"),
                                 header = TRUE,
                                 na.strings = c("NA", "NaN")))




colMeans(post_bias, na.rm = TRUE)
apply(post_bias, 2, sd, na.rm = TRUE)


post_len[is.infinite(post_len)] <- NA
colMeans(post_len, na.rm = TRUE)
apply(post_len, 2, sd, na.rm = TRUE)

colMeans(post_count, na.rm = TRUE)




1-colMeans(post_rej, na.rm = TRUE)



