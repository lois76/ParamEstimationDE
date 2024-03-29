//////////////////////////////////////////////////////////
///////////////     Experimental data      ///////////////
//////////////////////////////////////////////////////////

a = read("/scilab-scripts/Fig 1A_AIY Current-Clamp Trace.txt",-1,12);
A=a(2489:14988,2:$)*1000;
t=linspace(0,50,12500);
t0=0;
stim=[-15:5:35];

//Steady-state current
vecV=[-120:10:50];
Inf=[-13.1 -10.4 -7.92 -5.89 -4.11 -2.69 -1.02 0.0211 1.17 3.1 7.32 14.2 22.4 31.5 43.2 54.5 69.5 82.4];

//////////////////////////////////////////////////
///////////////    Cost function    //////////////
//////////////////////////////////////////////////

//Boltzmann function
function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

//Ca,t+Kir+K,p+L-model
function [Hdot]=HH11(t,x,pa)
    Hdot=zeros(4,1);
    Hdot(1)=(1/pa(22))*(-pa(1)*x(2)*x(3)*(x(1)-pa(5)) - pa(2)*xinf(x(1),pa(10),pa(14))*(x(1)-pa(6)) - pa(3)*x(4)*(x(1)-pa(6)) - pa(4)*(x(1)-pa(7)) + I)
    Hdot(2)=(xinf(x(1),pa(8),pa(12))-x(2))/pa(16)
    Hdot(3)=(xinf(x(1),pa(9),pa(13))-x(3))/pa(17)
    Hdot(4)=(xinf(x(1),pa(11),pa(15))-x(4))/pa(18)
endfunction

%ODEOPTIONS=[1,0,0,%inf,0,2,20000,12,5,0,-1,-1];

//Cost function voltage
function y=fct11(pa)
    c=0;
    condini = [-53; pa(19); pa(20); pa(21)]
    for i=1:11
        I=stim(i);
        x=ode(condini,t0,t,HH11); 
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
        e=e+(Inf(i)-(pa(1)*xinf(vecV(i),pa(8),pa(12))*xinf(vecV(i),pa(9),pa(13))*(vecV(i)-pa(5)) + pa(2)*xinf(vecV(i),pa(10),pa(14))*(vecV(i)-pa(6)) + pa(3)*xinf(vecV(i),pa(11),pa(15))*(vecV(i)-pa(6)) + pa(4)*(vecV(i)-pa(7))))^2
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

//////////////////////////////////////////////
/////////    Parameter estimation    /////////
//////////////////////////////////////////////

function [popInit, valInit, pop2500, val2500, pop5000, val5000, popFinal, valFinal]=simulation(NP,itermax,F,CR)
    
    D=22; 
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.1 0.1 0.1 0.1 20  -100 -90 -90 -90 -90 -90 1  -30 -30 1  0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001];
    Xmax=[50  50  50  50  150 -2   30  -2  -2  -2  -2  30 -1  -1  30 15     15     15     0.999 0.999 0.999 10];
    
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
    
    val=zeros(NP,2); // tableau avec le coût de chacun des individus. 1ère colonne = cout voltage. 2ème colonne = cout SS.
    
    for j=1:NP
        val(j,1)=fct11(pop(:,j))
        val(j,2)=W(pop(:,j))
    end
    
    // Save valInit
    valInit=val;

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
            tempvalVol = fct11(U(:,j));
            tempvalSS = W(U(:,j));
            if tempvalVol<tempval(j,1) & tempvalSS<=tempval(j,2) then
                tempPop(:,j) = U(:,j);
                tempval(j,1) = tempvalVol;
                tempval(j,2) = tempvalSS;
            end
            if (tempvalVol>tempval(j,1) & tempvalSS<=tempval(j,2)) | (tempvalVol<tempval(j,1) & tempvalSS>=tempval(j,2)) then
                tempPop=[tempPop U(:,j)]
                tempval=[tempval; [tempvalVol tempvalSS]]
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

        disp(iter);
        iter = iter + 1;
    end  //fin de la boucle while
    
    popFinal=pop;
    valFinal=val;
    disp(pop);
    disp(val);
endfunction

