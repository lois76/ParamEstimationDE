currentDate = getdate()
rand('seed',getdate('s')+currentDate(10))

exec("EstimationSteadyCurrent11.sce")

args=sciargs();
[valBest]=simulation(strtod(args(6)),strtod(args(7)),strtod(args(8)),strtod(args(9)));

csvWrite(valBest,"/Result/result.csv")
