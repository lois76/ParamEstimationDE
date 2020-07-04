//////////////////////////////////////////////////////////
///////////////     Experimental data      ///////////////
//////////////////////////////////////////////////////////

//Voltage
a = read("/home/loisse/Documents/FichierScilab/EstimationAFD/Fig1A_AFDCurrentClampTrace.txt",-1,12);
A=a(2489:14988,2:$)*1000;
t=linspace(0,50,12500);
t0=0;
stim=[-15:5:35];

//Steady-state current
vecV=[-110:10:50]
Inf=[-68.6 -49.5 -18.2 -5.06 2.19 3.37 2.52 2.68 5.97 14.6 33.4 60.2 85 114 152 208 254]

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

//Cost function voltage
function y=fct11(pa)
    c=0;
    condini = [-78; pa(19); pa(20); pa(21)]
    for i=1:11
        I=stim(i);
        x=ode(condini,t0,t,HH); 
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-A(k,i))*(V(k)-A(k,i))
        end
    end
    y=c/length(t);
endfunction

//Cost function steady-state current
function y=W(pa)
    e=0;
    for i=1:length(vecV)
        e=e+(Inf(i)-(pa(1).*xinf(vecV(i),pa(8),pa(12)).*(vecV(i)-pa(5)) + pa(2).*xinf(vecV(i),pa(9),pa(13)).*(vecV(i)-pa(6)) + pa(3).*xinf(vecV(i),pa(10),pa(14)).*xinf(vecV(i),pa(11),pa(15)).*(vecV(i)-pa(6)) + pa(4).*(vecV(i)-pa(7))))^2
    end
    y=e/length(vecV) 
endfunction

//////////////////////////////////////////////
/////////    Parameter estimation    /////////
//////////////////////////////////////////////

function [val, popFinal]=simulation(NP,itermax,F,CR)
    
    D=22; 
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.1 0.1 0.1 0.1 20  -100 -90 -90 -90 -90 -90 1  -30 1  -30 0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001];
    Xmax=[50  50  50  50  150 -2   30  -2  -2  -2  -2  30 -1  30 -1  15     15     15     0.999 0.999 0.999 10];
    
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
    
    val=zeros(NP,2); // tableau avec le coût de chacun des individus. 1ère colonne = cout voltage. 2ème colonne = cout SS.
    
    for j=1:NP
        val(j,1)=fct11(pop(:,j))
        val(j,2)=W(pop(:,j))
    end
    
    disp(val);
    
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
            if max(val)==min(val) then
                proba=ones(NP,1);
            else
                proba=exp(-5*(val-min(val))/(max(val)-min(val)));//donne un vecteur
            end
            proba=proba/sum(proba);

            r1=j; r2=j; r3=j;
            while (r1==r2 | r1==r3 | r2==r3 | r1==j | r2==j | r3==j)
                r1=find(1==grand(1,'mul',1,proba(1:$-1)))
                r2=find(1==grand(1,'mul',1,proba(1:$-1)))
                r3=find(1==grand(1,'mul',1,proba(1:$-1)))
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
            tempvalSS = W(U(:,j));
            if tempvalSS<=val(j,2) then
                tempvalVol = fct11(U(:,j));
                if tempvalVol<val(j,1) then
                    pop(:,j) = U(:,j);
                    val(j,1) = tempvalVol;
                    val(j,2) = tempvalSS;
                end
            end
        end
        
        disp(iter)
        iter = iter + 1;
    end  //fin de la boucle while

    // Détermination de l'indice du meilleur individu
//    bestIndex=1;
//    for b=2:NP
//        if val(b)<val(bestIndex) then bestIndex=b; end
//    end
//    valBest=val(bestIndex);
    
    // Sauvegarde du meilleur individu
//    bM = [];
//    bM = pop(:,bestIndex);
    popFinal=pop;
    disp(val);
endfunction
