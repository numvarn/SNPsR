library(SKAT)

# Config values
setwd("~/ResearchCode/SNPsR")

data(SKATBinary.example)
attach(SKATBinary.example)

data(SKAT.example)
attach(SKAT.example)


obj<-SKAT_Null_Model(y ~ x1 + x2, out_type="D")

SKATBinary(Z, obj, method.bin="ER.A")$p.value