args = sciargs();


vecV=[-110:10:50];
Inf=[-68.6 -49.5 -18.2 -5.06 2.19 3.37 2.52 2.68 5.97 14.6 33.4 60.2 85 114 152 208 254]

function y=xinf(V,V12,k)
    y=1 ./(1+exp((V12-V) ./k));
endfunction

// Fonction coût pour un modèle de type 2-1
function e=W21(pa)
    e=0;
    for i=1:length(vecV)
        e=e+(Inf(i)-(pa(1)*xinf(vecV(i),pa(7),pa(10))*(vecV(i)-pa(4))+pa(2)*xinf(vecV(i),pa(8),pa(11))*xinf(vecV(i),pa(9),pa(12))*(vecV(i)-pa(5))+pa(3)*(vecV(i)-pa(6))))^2
    end
endfunction

i = 6 
pa=[strtod(args(i)) strtod(args(i+1)) strtod(args(i+2)) strtod(args(i+3)) strtod(args(i+4)) strtod(args(i+5)) strtod(args(i+6)) strtod(args(i+7)) strtod(args(i+8)) strtod(args(i+9)) strtod(args(i+10)) strtod(args(i+11))];

[e]=W21(pa);


csvWrite(e,"/Result/result.csv")
