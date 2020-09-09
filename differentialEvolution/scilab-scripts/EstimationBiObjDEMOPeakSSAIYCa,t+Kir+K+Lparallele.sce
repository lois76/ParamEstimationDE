//////////////////////////////////////////////////////////
///////////////     Experimental data      ///////////////
//////////////////////////////////////////////////////////

t=linspace(0,50,12500);

//Steady-state current
vecV=[-120:10:50];
Inf=[-13.1 -10.4 -7.92 -5.89 -4.11 -2.69 -1.02 0.0211 1.17 3.1 7.32 14.2 22.4 31.5 43.2 54.5 69.5 82.4];
InfSD=[2.88 2.55 1.47 1.31 1.04 0.809 0.7 0.658 0.638 0.889 1.94 3.5 5.36 7.63 10.6 13.3 16 17.9]

//for i=1:length([-100:10:50])
//    plot(vecV(i),Inf(i),"r.")
//    plot(vecV(i),Inf(i)+InfSD(i),"r+")
//    plot(vecV(i),Inf(i)-InfSD(i),"r+")
//end

//Peak current
vecVpk=[-120:10:50]
pk=[-11.9 -9.42 -7.39 -4.88 -3.16 -1.48 0.266 1.67 2.56 4.19 8.25 15.7 24.5 36.4 47.8 62.2 78.2 91.9]
pkSD=[1.7 1.46 0.979 0.834 0.807 0.567 0.559 0.751 1.06 1.14 1.89 3.29 5.26 8.09 11.5 14.6 18.1 21.3]

//for i=1:length([-100:10:50])
//    plot(vecVpk(i),pk(i),"r.")
//    plot(vecVpk(i),pk(i)+InfSD(i),"r+")
//    plot(vecVpk(i),pk(i)-InfSD(i),"r+")
//end

///////////////////////////////////////////////////////////////////////
///////////////    Steady-state current cost function    //////////////
///////////////////////////////////////////////////////////////////////

//Boltzmann function
function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

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



