columnnum=2
ASVTADMMErr=ASVTADMMErrresultsCollection[:,columnnum]
ASVTErr=ASVTErrresultsCollection[:,columnnum]
lm=(length(ASVTErr)-1)
using Plots
x = 1:lm;
y=zeros(lm,4)
y[:,1]= ASVTErr[1:lm];
y[:,2]= ASVTErr[lm+1]*ones(lm,1)
y[:,3]= ASVTADMMErr[1:lm];
y[:,4]= ASVTADMMErr[lm+1]*ones(lm,1)
Plots.plot(x,y)
