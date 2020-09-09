//////////////////////////////////////////////////////////
///////////////     Experimental data      ///////////////
//////////////////////////////////////////////////////////

t=linspace(0,50,12500);
//tpeak=t(1:251);

//Steady-state current
vecV=[-110:10:50]
Inf=[-68.6 -49.5 -18.2 -5.06 2.19 3.37 2.52 2.68 5.97 14.6 33.4 60.2 85 114 152 208 254]
InfSD=[1 8.65 0.636 1.31 1.83 1.46 0.814 0.455 0.613 2.63 7.71 14.7 22.3 27.4 44.1 73.7 97.6]

//for i=1:length([-100:10:50])
//    plot(vecV(i),Inf(i),"r.")
//    plot(vecV(i),Inf(i)+InfSD(i),"r+")
//    plot(vecV(i),Inf(i)-InfSD(i),"r+")
//end

//Peak current
vecVpk=[-120:10:50]
pk=[-75.2 -55.7 -34.6 -16.3 -2.11 3.91 5.36 6.15 6.58 9.71 30.2 65.7 109 159 213 260 312 359]
pkSD=[3.15 1.45 0.674 1.08 1.44 1.94 2.68 2.03 1.5 1.22 4.11 10.7 18.1 23.1 31.7 37.2 43.7 50]

//for i=1:length([-100:10:50])
//    plot(vecVpk(i),pk(i),"r.")
//    plot(vecVpk(i),pk(i)+InfSD(i),"r+")
//    plot(vecVpk(i),pk(i)-InfSD(i),"r+")
//end

//////////////////////////////////////////////////////////////
///////////////    Cost function peak current   //////////////
//////////////////////////////////////////////////////////////

//Boltzmann function
function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

//Solution x(t) where x=m,h
function y=x(VH,V12,k,tx,x0,t)
    y=xinf(VH,V12,k)+(x0-xinf(VH,V12,k))*exp(-t/tx)
endfunction

//Current equation
function y=Iest(VH,pa,t)
    y = pa(1).*x(VH,pa(8),pa(12),pa(16),pa(19),t).*(VH-pa(5)) + pa(2)*xinf(VH,pa(9),pa(13))*(VH-pa(6)) + pa(3).*x(VH,pa(10),pa(14),pa(17),pa(20),t).*x(VH,pa(11),pa(15),pa(18),pa(21),t).*(VH-pa(6)) + pa(4).*(VH-pa(7))
endfunction

//Cost function peak current
function y=Wpeak(pa)
    e=0;
    for i=1:length(vecVpk)
        Ipeak=0;
        tmp=0;
        Ipeak=max(Iest(vecVpk(i),pa,t(1:251)));
        tmp=(pk(i)-Ipeak)^2;
        tmp=tmp/pkSD(i);
        e=e+tmp;
    end
    y=e/length(vecVpk) 
endfunction

//////////////////////////////////////////////////////////////////////
///////////////    Cost function steady-state current   //////////////
//////////////////////////////////////////////////////////////////////

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
/////////    Crowding Sorting and Domination Front    /////////
///////////////////////////////////////////////////////////////

function [Front]=NDS(A)
    dominationCount = zeros(size(A,'r'),1)
    S=list(); // S(1) contiendra les indices des solutions que la solution i domine
    Front=list(); // F(1) contiendra les indices des solutions du front 1, F(2) les indices des solutions des fronts 2, etc...
    for i=1:size(A,'r')
        Stmp=[]; // vecteur vide qui va contenir les solutions que la solution i domine
        for j=1:size(A,'r')
            if i~=j then 
                // nombre de solutions qui domine la solution i
                if A(i,1)>A(j,1) & A(i,2)>A(j,2) then
                    dominationCount(i) = dominationCount(i) + 1;
                end
                // Ensemble de solution que la solution i domine
                if A(i,1)<A(j,1) & A(i,2)<A(j,2) then
                    Stmp=[Stmp j]
                end
            end
        end
        S(i)=Stmp
    end
    Front(1)=find(0==dominationCount); // indice des solutions faisant partie du best front
    
    // Il faudra "set a front counter m=1 pour itérer sur tous les fronts pour des cas plus complexes et un while.. for i=
    m=1
    while Front(m)~=[] // tant que le front m est non vide alors...
        Q=[];
        for i=Front(m)
            for j=S(i)
                dominationCount(j) = dominationCount(j) - 1;
                if dominationCount(j)==0 then
                    Q=[Q j];
                end
            end
        end
        m=m+1;
        Front(m)=Q;
    end
endfunction


function [d]=crowdingSorting(A)
    l=size(A, 1) // l=nombre d'individus dans A (qui est l'ensemble des fonction coûts du dernier front)
    d = zeros(l, 1);
    for m=1:2 // pour chaque fonction coûts m. Ici m=2 car seulement 2 fonctions coûts
        [tmp, Index] = gsort(A(:, m)); // Step C2 : sort the set in ascendant order of magnitude
//        ////pause;
        d(Index(1)) = %inf;
        d(Index(l)) = %inf;
        fmax = max(A(:, m));
        fmin = min(A(:, m));
//        ////pause;
        for j=2:l-1
            d(Index(j)) = d(Index(j)) + abs(tmp(j+1) - tmp(j-1)) / (fmax - fmin);
        end
//        //pause;
    end
endfunction


///////////////////////////////////////////////////
/////////    Estimation des paramètres    /////////
///////////////////////////////////////////////////

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
    // Save popInit
    popInit=pop;
    
//    //pause;
    ///////////////////////////////////////
    //// Initial population evaluation ////
    ///////////////////////////////////////
    
    val=zeros(NP,2); // tableau avec le coût de chacun des individus. 1ère colonne = cout SS. 2ème colonne = cout peak.
    
    for j=1:NP
        val(j,1)=WSS(pop(:,j))
        val(j,2)=Wpeak(pop(:,j))
    end
    
    // Save valInit
    valInit=val;
    disp(valInit);
//    //pause;
    ////////////////////////
    //// Étape suivante ////
    ////////////////////////
     
    iter=1; // nombre d'itération
    U=zeros(D,NP); // Vecteur intermédiaire perturbé (mutation + crossover)
    tempvalVol=0;
    tempvalSS=0;
    while iter<itermax
        for j=1:NP
            // ======= Construction de la matrice U = variation différentielle + crossover =======

            // ========= Tirage aléatoire de 3 entiers distincts r1, r2 et r3 et différents de j ========
            r1=j; r2=j; r3=j;//////////////////////////////////////
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
        // Ajout des enfants U dans la pop si ils dominent les parents ou si ils sont non-dominés. Donc |pop|>NP
        tempPop=pop; // tempPop sera la population modifié et augmenté avec les enfants qui dominent les parents ou enfant+parent qui sont non dominés        
        tempval=val;
        
        for j=1:NP
            tempvalSS = WSS(U(:,j));
            tempvalpeak = Wpeak(U(:,j));
            if tempvalSS<tempval(j,1) & tempvalpeak<=tempval(j,2) then
                tempPop(:,j) = U(:,j);
                tempval(j,1) = tempvalSS;
                tempval(j,2) = tempvalpeak;
            end
            if (tempvalSS>tempval(j,1) & tempvalpeak<=tempval(j,2)) | (tempvalSS<tempval(j,1) & tempvalpeak>=tempval(j,2)) then
                tempPop=[tempPop U(:,j)]
                tempval=[tempval; [tempvalSS tempvalpeak]]
            end
        end
        
        // Front ranking de tempPop > NP
        [Front]=NDS(tempval);
        
        // Intégration des fronts possibles dans la pop
        pop=[];
        val=[];
        k=1;
        while (size(pop,2)+length(Front(k)))<NP
            for i=1:length(Front(k))
                pop=[pop tempPop(:,Front(k)(i))];
                val=[val; tempval(Front(k)(i),:)];
            end
            k=k+1;
        end
       
        // Calcul de la distance de crowding du dernier front F(k) qui doit être tronqué
        lastFront=[];
        for i=1:length(Front(k))
            lastFront=[lastFront; tempval(Front(k)(i),:)];
        end
        
        cs=crowdingSorting(lastFront);//Asignation d'une distance de crowding
        
        // Intégration des individus du dernier front selon leur distance de crowding
        [osef, indice]=gsort(cs);
        
        n=1;
        while size(pop,2)<NP
            pop=[pop tempPop(:,Front(k)(indice(n)))];
            val=[val; tempval(Front(k)(indice(n)),:)];
            n=n+1;
        end
        
        if iter==8000 then
            disp(pop);
            disp(val);
            pop2500=pop;
            val2500=val;
        end
        if iter==17000 then
            disp(pop);
            disp(val);
            pop5000=pop;
            val5000=val;
        end

        disp(iter);
        iter = iter + 1;
    end  //fin de la boucle while
    
    popFinal=pop;
    valFinal=val;
    disp(pop);
    disp(val);
endfunction



