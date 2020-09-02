//////////////////////////////////////////////////////////
///////////////     Experimental data      ///////////////
//////////////////////////////////////////////////////////

a = read("/home/naudin/Documents/Fig1A_AFDCurrentClampTrace.txt",-1,12);
//a = read("/scilab-scripts/Fig 1A_AFD Current-Clamp Trace.txt",-1,12);
A=a(2489:14988,2:$)*1000;
t=linspace(0,50,12500);
t0=0;
stim=[-15:5:35];

//Steady-state current
vecV=[-110:10:50]
Inf=[-68.6 -49.5 -18.2 -5.06 2.19 3.37 2.52 2.68 5.97 14.6 33.4 60.2 85 114 152 208 254]
InfSD=[1 8.65 0.636 1.31 1.83 1.46 0.814 0.455 0.613 2.63 7.71 14.7 22.3 27.4 44.1 73.7 97.6]

//Peak current
vecVpk=[-120:10:50]
pk=[-75.2 -55.7 -34.6 -16.3 -2.11 3.91 5.36 6.15 6.58 9.71 30.2 65.7 109 159 213 260 312 359]
pkSD=[3.15 1.45 0.674 1.08 1.44 1.94 2.68 2.03 1.5 1.22 4.11 10.7 18.1 23.1 31.7 37.2 43.7 50]

//////////////////////////////////////////////////
///////////////    Cost function    //////////////
//////////////////////////////////////////////////

//Boltzmann function
function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

//Ca,p+Kir+K,t+L-model
function [Hdot]=HH(t,x,pa)
    Hdot=zeros(4,1);
    Hdot(1)=(1/pa(22))*(-pa(1)*x(2)*(x(1)-pa(5)) - pa(2)*xinf(x(1),pa(9),pa(13))*(x(1)-pa(6)) - pa(3)*x(3)*x(4)*(x(1)-pa(6)) - pa(4)*(x(1)-pa(7)) + I)
    Hdot(2)=(xinf(x(1),pa(8),pa(12))-x(2))/pa(16)
    Hdot(3)=(xinf(x(1),pa(10),pa(14))-x(3))/pa(17)
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

// Définition de l'écart type sur V pour chaque I
dev=[]
for i=1:11
    dev=[dev sigma(A(7000:$,i))]
end

%ODEOPTIONS=[1,0,0,%inf,0,2,20000,12,5,0,-1,-1];

//Cost function voltage
function y=fct11(pa)
    tmp=0;
    condini = [-78; pa(19); pa(20); pa(21)]
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
        tmp=(Inf(i)-(pa(1).*xinf(vecV(i),pa(8),pa(12)).*(vecV(i)-pa(5)) + pa(2).*xinf(vecV(i),pa(9),pa(13)).*(vecV(i)-pa(6)) + pa(3).*xinf(vecV(i),pa(10),pa(14)).*xinf(vecV(i),pa(11),pa(15)).*(vecV(i)-pa(6)) + pa(4).*(vecV(i)-pa(7))))^2
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
    y = pa(1).*x(VH,pa(8),pa(12),pa(16),pa(19),t).*(VH-pa(5)) + pa(2)*xinf(VH,pa(9),pa(13))*(VH-pa(6)) + pa(3).*x(VH,pa(10),pa(14),pa(17),pa(20),t).*x(VH,pa(11),pa(15),pa(18),pa(21),t).*(VH-pa(6)) + pa(4).*(VH-pa(7))
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

    Xmin=[0.0001 0.0001 0.0001 0.0001 20  -100 -90 -90 -90 -90 -90 1  -30 1  -30 0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001];
    Xmax=[50     50     50     50     150 -2   30  -2  -2  -2  -2  30 -1  30 -1  20     20     20     0.999 0.999 0.999 10];
    
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
