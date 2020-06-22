currentDate = getdate()
rand('seed',getdate('s')+currentDate(10))

args=sciargs();

exec("/scilab-scripts/"+args(6))

[bM, valBest]=simulation(strtod(args(7)),strtod(args(8)),strtod(args(9)),strtod(args(10)));

csvWrite(valBest,"/Result/result.csv")
csvWrite(bM,"/Result/bM.csv")
csvWrite(costVec,"/Result/costVec.csv")

