library(NBZIMM)

data(Romero)
names(Romero)

otu = Romero$OTU; dim(otu)
sam = Romero$SampleData; dim(sam)
colnames(sam)

N = sam[, "Total.Read.Counts"]        
Days = sam$GA_Days; Days = scale(Days)
Age = sam$Age; Age = scale(Age)
Race = sam$Race
preg = sam$pregnant; table(preg)
subject = sam[, "Subect_ID"]; table(subject)

#y = otu[, 1]

#f1 = glmm.nb(y ~ Days + Age + Race + preg + offset(log(N)), random = ~ 1 | subject)
#summary(f1)
#fixed(f1)

#For all taxa, we can analyze them using the fuction mms for all four above models:

# The first model:
f2 = mms(y = Romero$OTU, fixed = ~  Days + Age + Race + preg + offset(log(N)), 
        random = ~ 1 | subject, data = Romero$SampleData,
        min.p = 0.2, method = "nb")

fixed(f2)
       
