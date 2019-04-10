library(vegan)
m1 <- cca(beta ~ ., data = mut_sig)
m0 <- cca(beta ~ 1, data = mut_sig)
m <- step(m0, scope=formula(m1), test="perm")
mback <- step(m1, test="perm")


mval = mval + abs(min(mval))
m1 <- cca(mval ~ ., data = mut_sig)
m0 <- cca(mval ~ 1, data = mut_sig)
m <- step(m0, scope=formula(m1), test="perm")
mback <- step(m1, test="perm")
