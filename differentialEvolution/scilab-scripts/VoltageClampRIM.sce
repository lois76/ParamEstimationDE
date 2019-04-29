args = sciargs();


// Récupération des données expérimentales sous forme de matrice
a = read("/scilab-scripts/DataRIMVoltageClamp.txt",-1,17);
A=a(239:1488,2:$)*1E12;
t=linspace(0,50,1250);

////////////////////////////////////////////////////////////////
//// Estimation de [tx1 tx2 tx3 tx4 x1(0) x2(0) x3(0) x4(0) ////
////////////////////////////////////////////////////////////////

gCa=0.175; gK=30; gL=0.366;
ECa=113.32; EK=-38.646; EL=-68.808;
V12x1=-34.411; V12x2=-37.415; V12x3=-69.735;
kx1=16.995; kx2=5.225; kx3=-5.658;

function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

function y=x(VH,V12,k,tx,x0,t)
    y=xinf(VH,V12,k)+(x0-xinf(VH,V12,k))*exp(-t/tx)
endfunction

function y=Iest(VH,pa,t)
    y=gCa.*x(VH,V12x1,kx1,pa(1),pa(4),t).*(VH-ECa) + gK.*x(VH,V12x2,kx2,pa(2),pa(5),t).*x(VH,V12x3,kx3,pa(3),pa(6),t).*(VH-EK) + gL.*(VH-EL)
endfunction


// Fonction coût
VH=[-100:10:50];

function y=W(pa)
    c=0;
    for i=1:length(VH)
        for k=1:length(t)
            c=c+(Iest(VH(i),pa,t(k))-A(k,i)).^2
        end
    end
    y=c/length(t);
endfunction

i = 6;
pa=[strtod(args(i)) strtod(args(i+1)) strtod(args(i+2)) strtod(args(i+3)) strtod(args(i+4)) strtod(args(i+5))];

[e]=W(pa);


csvWrite(e,"/Result/result.csv")


