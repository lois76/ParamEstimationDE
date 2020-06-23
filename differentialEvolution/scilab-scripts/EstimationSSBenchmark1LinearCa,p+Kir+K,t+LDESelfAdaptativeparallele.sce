////////////////////////////////////////////////////////////////////////
///////////////     Benchmark Linear Ca,p+Kir+K,t+L      ///////////////
////////////////////////////////////////////////////////////////////////

// Let bc be the vector storing the parameters of benchmark with monotonic steady-state current for a ICa,p+IKir+IK,t+IL-model. As you can see below, gCa is the first parameter of bc, gKir the second, gK the third, etc...  

bc=[0.3 
    0.12  
    0.86  
    0.44  
    35.3  
  -77.1  
  -64.23  
  -29.75  
  -74.6   
  -48.9  
  -79.6  
    29.1  
  -15.1  
    28.71  
  -21.0]

gCa=bc(1); gKir=bc(2); gK=bc(3); gL=bc(4);
ECa=bc(5); EK=bc(6); EL=bc(7);
V12mCa=bc(8); V12hKir=bc(9); V12mK=bc(10); V12hK=bc(11);
kmCa=bc(12); kKir=bc(13); kmK=bc(14); khK=bc(15);

/////////////////////////////////////////////////////////////////
///////////////    Plot of steady-state current    //////////////
/////////////////////////////////////////////////////////////////

vecV=[-100:10:50]

function y=xinf(V,V12,k)
    y=1 ./(1+exp((V12-V) ./k));
endfunction

function y=Iinf(VH)
    y = gCa.*xinf(VH,V12mCa,kmCa).*(VH-ECa) + gKir.*xinf(VH,V12hKir,kKir).*(VH-EK) + gK.*xinf(VH,V12mK,kmK).*xinf(VH,V12hK,khK).*(VH-EK) + gL.*(VH-EL)
endfunction

/////////////////////////////////////////////////////////////
///////////////    Fonction coût algorithme    //////////////
/////////////////////////////////////////////////////////////


function y=W(pa)
    e=0;
    for i=1:length(vecV)
        e=e+(Iinf(vecV(i))-(pa(1).*xinf(vecV(i),pa(8),pa(12)).*(vecV(i)-pa(5)) + pa(2).*xinf(vecV(i),pa(9),pa(13)).*(vecV(i)-pa(6)) + pa(3).*xinf(vecV(i),pa(10),pa(14)).*xinf(vecV(i),pa(11),pa(15)).*(vecV(i)-pa(6)) + pa(4).*(vecV(i)-pa(7))))^2
    end
    y=e/length(vecV) 
endfunction


//////////////////////////////////////////////
/////////    Parameter estimation    /////////
//////////////////////////////////////////////

Fl=0.1; Fu=0.9;
t1=0.1; t2=0.1;

function [bM, valBest, costVec]=simulation(NP,itermax,F,CR)
    
    D=17;//15+2 because of F and CR parameters included in individuals
    costVec=zeros(1,itermax);
    pop=zeros(D,NP);
    
    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.1 0.1 0.1 0.1 20  -100 -90 -90 -90 -90 -90 1  -30 1  -30 0.1 0];
    Xmax=[50  50  50  50  150 -2   30  -2  -2  -2  -2  30 -1  30 -1  1   1];
    
    /////////////////////////////////////////
    //// Initialisation de ma population ////
    /////////////////////////////////////////

    for j=1:NP
        for i=1:D
            pop(i,j)=Xmin(i)+(Xmax(i)-Xmin(i))*rand();
        end
    end
    
    ///////////////////////////////////////////////////////////////////
    //// Évaluation du meilleur individu de ma première génération ////
    ///////////////////////////////////////////////////////////////////
    
    val=zeros(NP,1); // tableau avec le coût de chacun des individus
    
    for j=1:NP
        val(j)=W(pop(:,j))
    end
    
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
    costVec(1)=val(bestIndex);
    
    vecF=zeros(1,NP)
    vecCR=zeros(1,NP)
    vecF(1)=pop(D-1,bestIndex)
    vecCR(1)=pop(D,bestIndex)
    
    ////////////////////////
    //// Étape suivante ////
    ////////////////////////
     
    iter=1; // nombre d'itération
    
    U=zeros(D,NP); // Matrice intermédiaire perturbé (mutation + crossover)
    for j=1:NP
        U(D-1,j)=pop(D-1,j)
        U(D,j)=pop(D,j)
    end

    tempval=0;
    
    while iter<itermax
        for j=1:NP
            // ======= Construction de la matrice U = variation différentielle + crossover =======

            // ======== Differential variation param DE for new parent vector =======
            if rand()<t1 then
                U(D-1,j)=Fl+rand()*Fu
            end

            if rand()<t2 then
                U(D,j)=rand();
            end
            
            // ========= Tirage aléatoire de 3 entiers distincts r1, r2 et r3 et différents de j ========
            // ======== Variation différentielle param modèle =======
            r1=j; r2=j; r3=j;
            while (r1==r2 | r1==r3 | r2==r3 | r1==j | r2==j | r3==j)
                r1=floor(1+NP*rand());
                r2=floor(1+NP*rand());
                r3=floor(1+NP*rand());
            end

            V(1:D-2)=pop(1:D-2,r1) + U(D-1,j)*(pop(1:D-2,r2)-pop(1:D-2,r3));

            // ======== Constraints ========
            for i=1:D-2
                if V(i)<=Xmin(i) then V(i)=Xmin(i);
                elseif V(i)>Xmax(i) then V(i)=Xmax(i);
                end
            end
            
            // ======== Crossover ========
            for i=1:D-2
                if rand()<U(D,j) then
                    U(i,j)=V(i);
                else
                    U(i,j)=pop(i,j);
                end
            end
        end // fin for j=1:NP
    
    // ======== Sélection ========
        for j=1:NP
            tempval=W(U(:,j));

            if tempval<=val(j) then
                pop(:,j) = U(:,j);
                val(j) = tempval;
            end
        end
        disp(iter)
        iter = iter + 1;
        bestIndex=1;
        for b=2:NP
            if val(b)<val(bestIndex) then bestIndex=b; end
        end
        costVec(iter)=val(bestIndex);
        vecF(iter)=pop(D-1,bestIndex);
        vecCR(iter)=pop(D,bestIndex);
    end  //fin de la boucle while
    
    // Détermination de l'indice du meilleur individu
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
//    disp(bestIndex);
    
    // Sauvegarde du meilleur individu
    bM = [];
    bM = pop(:,bestIndex);
    
    valBest=val(bestIndex);
    
    disp(val);
//    disp(bM);
    disp(val(bestIndex));
    
//    iterVec=1:1:itermax;
//    plot(iterVec,costVec,2)
endfunction


