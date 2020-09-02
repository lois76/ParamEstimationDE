//////////////////////////////////////////////////////////
///////////////     Experimental data      ///////////////
//////////////////////////////////////////////////////////

//a = read("/home/naudin/Documents/FichierScilab/Fourre tout/Fig1A_AIYCurrentClampTrace2.txt",-1,12);
a = read("/scilab-scripts/Fig 1A_AIY Current-Clamp Trace.txt",-1,12);
A=a(2489:14988,2:$)*1000;
t=linspace(0,50,12500);
t0=0;
stim=[-15:5:35];

//Steady-state current
vecV=[-120:10:50];
Inf=[-13.1 -10.4 -7.92 -5.89 -4.11 -2.69 -1.02 0.0211 1.17 3.1 7.32 14.2 22.4 31.5 43.2 54.5 69.5 82.4];
InfSD=[2.88 2.55 1.47 1.31 1.04 0.809 0.7 0.658 0.638 0.889 1.94 3.5 5.36 7.63 10.6 13.3 16 17.9]

//Peak current
vecVpk=[-120:10:50]
pk=[-11.9 -9.42 -7.39 -4.88 -3.16 -1.48 0.266 1.67 2.56 4.19 8.25 15.7 24.5 36.4 47.8 62.2 78.2 91.9]
pkSD=[1.7 1.46 0.979 0.834 0.807 0.567 0.559 0.751 1.06 1.14 1.89 3.29 5.26 8.09 11.5 14.6 18.1 21.3]

//////////////////////////////////////////////////
///////////////    Cost function    //////////////
//////////////////////////////////////////////////

//Boltzmann function
function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

//Ca,t+Kir+K,p+L-model
function [Hdot]=HH(t,x,pa)
    Hdot=zeros(4,1);
    Hdot(1)=(1/pa(22))*(-pa(1)*x(2)*x(3)*(x(1)-pa(5)) - pa(2)*xinf(x(1),pa(10),pa(14))*(x(1)-pa(6)) - pa(3)*x(4)*(x(1)-pa(6)) - pa(4)*(x(1)-pa(7)) + I)
    Hdot(2)=(xinf(x(1),pa(8),pa(12))-x(2))/pa(16)
    Hdot(3)=(xinf(x(1),pa(9),pa(13))-x(3))/pa(17)
    Hdot(4)=(xinf(x(1),pa(11),pa(15))-x(4))/pa(18)
endfunction

//Fonction qui calcule l'écart type (standard deviation)
function y=sigma(v)
    s=0;
    moy=mean(v);
    for i=1:length(v)
        s=s+(v(i)-moy)^2
    end
    y=sqrt(s/(length(v)-1));
endfunction

//Noise level (standard deviation) for each I
dev=[]
dev1=sigma(A(500:$,1));
dev=[dev1]
for i=2:11
    dev=[dev sigma(A(5000:$,i))]
end

%ODEOPTIONS=[1,0,0,%inf,0,2,20000,12,5,0,-1,-1];

//Cost function voltage
function y=fct11(pa)
    tmp=0;
    condini = [-53; pa(19); pa(20); pa(21)]
    for i=1:11
        c=0;
        I=stim(i);
        x=ode(condini,t0,t,HH); 
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-A(k,i))*(V(k)-A(k,i));
        end
        c=sqrt(c/length(t))/dev(i);
        tmp=tmp+c;
    end
    y=tmp/11;
endfunction

//Cost function steady-state current
///////////////////////////////////////////////////////////////////////
///////////////    Steady-state current cost function    //////////////
///////////////////////////////////////////////////////////////////////

function y=WSS(pa)
    e=0;
    for i=1:length(vecV)
        tmp=0;
        tmp=(Inf(i)-(pa(1)*xinf(vecV(i),pa(8),pa(12))*xinf(vecV(i),pa(9),pa(13))*(vecV(i)-pa(5)) + pa(2)*xinf(vecV(i),pa(10),pa(14))*(vecV(i)-pa(6)) + pa(3)*xinf(vecV(i),pa(11),pa(15))*(vecV(i)-pa(6)) + pa(4)*(vecV(i)-pa(7))))^2
        tmp=tmp/InfSD(i)
        e=e+tmp;
    end
    y=e/length(vecV) 
endfunction

///////////////////////////////////////////////////////////////
///////////////    Peak current cost function    //////////////
///////////////////////////////////////////////////////////////

//Solution x(t) where x=m,h
function y=x(VH,V12,k,tx,x0,t)
    y=xinf(VH,V12,k)+(x0-xinf(VH,V12,k))*exp(-t/tx)
endfunction

//Current equation
function y=Iest(VH,pa,t)
    y = pa(1).*x(VH,pa(8),pa(12),pa(16),pa(19),t).*x(VH,pa(9),pa(13),pa(17),pa(20),t).*(VH-pa(5)) + pa(2)*xinf(VH,pa(10),pa(14))*(VH-pa(6)) + pa(3).*x(VH,pa(11),pa(15),pa(18),pa(21),t).*(VH-pa(6)) + pa(4).*(VH-pa(7))
endfunction

//Cost function peak current
function y=Wpk(pa)
    e=0;
    for i=1:length(vecVpk)
        Ipeak=[];
        tmp=0;
        Ipeak=max(Iest(vecVpk(i),pa,t(1:251)));
        tmp=(pk(i)-Ipeak)^2;
        tmp=tmp/pkSD(i);
        e=e+tmp;
    end
    y=e/length(vecVpk) 
endfunction


//////////////////////////////////////////////
/////////    Parameter estimation    /////////
//////////////////////////////////////////////

function [popInit, valInit, pop2500, val2500, pop5000, val5000, popFinal, valFinal]=simulation(NP,itermax,F,CR)
    
    D=22; 
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.0001 0.0001 0.0001 0.0001 20  -100 -90 -90 -90 -90 -90 1  -30 -30 1  0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001];
    Xmax=[50     50     50     50     150 -2   30  -2  -2  -2  -2  30 -1  -1  30 20     20     20     0.999 0.999 0.999 10];
    
    ////////////////////////////////////
    ////  Population initialization ////
    ////////////////////////////////////

    for j=1:NP
        for i=1:D
            pop(i,j)=Xmin(i)+(Xmax(i)-Xmin(i))*rand();
        end
    end
    
    ///////////////////////////////////////
    //// Initial population evaluation ////
    ///////////////////////////////////////
    
    val=zeros(NP,3); // Array with cost function of each individual. 1st column = voltage cost; 2nd column = steady-state cost; 3st column = peak cost
    
    for j=1:NP
        val(j,1)=fct11(pop(:,j));
        val(j,2)=WSS(pop(:,j));
        val(j,3)=Wpk(pop(:,j));
    end
    
    // Save valInit
    valInit=val;
    disp(valInit);
    
    ////////////////////////
    //// Étape suivante ////
    ////////////////////////
     
    iter=1; // nombre d'itération
    U=zeros(D,NP); // Vecteur intermédiaire perturbé (mutation + crossover)
    tempval=0;
    while iter<itermax
        for j=1:NP
            // ======= Construction de la matrice U = variation différentielle + crossover =======

            // ========= Tirage aléatoire de 3 entiers distincts r1, r2 et r3 et différents de j ========
            r1=j; r2=j; r3=j;
            while (r1==r2 | r1==r3 | r2==r3 | r1==j | r2==j | r3==j)
                r1=floor(1+NP*rand());
                r2=floor(1+NP*rand());
                r3=floor(1+NP*rand());
            end
            // ======== Variation différentielle =======
            V=pop(:,r1) + F*(pop(:,r2)-pop(:,r3));
            
            // ======== Contraintes ========
            for i=1:length(Xmin)
                if V(i)<=Xmin(i) then V(i)=Xmin(i);
                elseif V(i)>Xmax(i) then V(i)=Xmax(i);
                end
            end
            // ======== Crossover ========
            for i=1:D
                if rand()<CR then
                    U(i,j)=V(i);
                else
                    U(i,j)=pop(i,j);
                end
            end
        end // fin for j=1:NP

        // ======== Sélection ========
        for j=1:NP
            tempvalSS = WSS(U(:,j));
            if tempvalSS <= val(j,2) then
                tempvalPeak = Wpk(U(:,j));
                if tempvalPeak<val(j,3) then
                    tempvalVol = fct11(U(:,j));
                    if tempvalVol<val(j,1) then
                        pop(:,j) = U(:,j);
                        val(j,1) = tempvalVol;
                        val(j,2) = tempvalSS;
                        val(j,3) = tempvalPeak;
                    end
                end
            end
        end

        if iter==2500 then
            disp(pop);
            disp(val);
            pop2500=pop;
            val2500=val;
        end
        if iter==5000 then
            disp(pop);
            disp(val);
            pop5000=pop;
            val5000=val;
        end

        if (iter==3 | iter==300 | iter==600 | iter==900 | iter==1200 | iter==1500 | iter==2000 | iter==3000 | iter==3500 | iter==4000 | iter==4500 | iter==5500) then
            disp(pop);
            disp(val);
        end
        
        disp(iter)
        iter = iter + 1;
    end  //fin de la boucle while

    popFinal=pop;
    valFinal=val;
    disp(pop);
    disp(val);
endfunction

[popInit, valInit, pop2500, val2500, pop5000, val5000, popFinal, valFinal]=simulation(40,300,0.5,0.9)
