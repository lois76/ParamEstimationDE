////////////////////////////////////////////////////////////////////////
///////////////     Benchmark Linear Ca,p+Kir+K,t+L      ///////////////
////////////////////////////////////////////////////////////////////////

// Let bc be the vector storing the parameters of benchmark with monotonic steady-state current for a ICa,p+IKir+IK,t+IL-model. As you can see below, gCa is the first parameter of bc, gKir the second, gK the third, etc...  
  
bc=[4.7  
    1.3  
    2.6  
    1.6  
    59.7 
  -98.9  
  -88.4  
  -3.4  
  -72.6  
  -8.2  
  -55.4  
    28.6  
  -19.6 
    14.1  
  -15.4]

gCa=bc(1); gKir=bc(2); gK=bc(3); gL=bc(4);
ECa=bc(5); EK=bc(6); EL=bc(7);
V12mCa=bc(8); V12hKir=bc(9); V12mK=bc(10); V12hK=bc(11);
kmCa=bc(12); kKir=bc(13); kmK=bc(14); khK=bc(15);


/////////////////////////////////////////////////////////////////
///////////////    Plot of steady-state current    //////////////
/////////////////////////////////////////////////////////////////

vecV=[-110:10:50]

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

function [bM, valBest]=simulation(NP,itermax,F,CR)
    
    D=15;
//    costVec=zeros(1,itermax);
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.1 0.1 0.1 0.1 20  -100 -90 -90 -90 -90 -90 1  -30 1  -30];
    Xmax=[50  50  50  50  150 -2   30  -2  -2  -2  -2  30 -1  30 -1];
    
    /////////////////////////////////////////
    //// Initialisation de ma population ////
    /////////////////////////////////////////

    for j=1:NP
        for i=1:D
            pop(i,j)=Xmin(i)+(Xmax(i)-Xmin(i))*rand();
        end
    end
    

    //////////////////////////////////////////////////////////////
    //// Évaluation du meilleur individu après initialisation ////
    //////////////////////////////////////////////////////////////
    
    val=zeros(NP,1); // tableau avec le coût de chacun des individus
    
    for j=1:NP
        val(j)=W(pop(:,j))
    end
    
    disp(val);
    
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
    costVec(1)=val(bestIndex);
    
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
            r1=j; r2=j; r3=j;//////////////////////////////////////
            while (r1==r2 | r1==r3 | r2==r3 | r1==j | r2==j | r3==j)
                r1=floor(1+NP*rand());
                r2=floor(1+NP*rand());
                r3=floor(1+NP*rand());
            end
            // ======== Variation différentielle =======
            V=pop(:,r1) + F*(pop(:,r2)-pop(:,r3));
            
            if V(1)<=Xmin(1) then V(1)=Xmin(1);
            elseif V(1)>Xmax(1) then V(1)=Xmax(1);
            end
            if V(2)<=Xmin(2) then V(2)=Xmin(2);
            elseif V(2)>Xmax(2) then V(2)=Xmax(2);
            end
            if V(3)<=Xmin(3) then V(3)=Xmin(3);
            elseif V(3)>Xmax(3) then V(3)=Xmax(3);
            end
            if V(4)<=Xmin(4) then V(4)=Xmin(4);
            elseif V(4)>Xmax(4) then V(4)=Xmax(4);
            end
            if V(5)<=Xmin(5) then V(5)=Xmin(5);
            elseif V(5)>Xmax(5) then V(5)=Xmax(5);
            end
            if V(6)<=Xmin(6) then V(6)=Xmin(6);
            elseif V(6)>Xmax(6) then V(6)=Xmax(6);
            end
            if V(7)<=Xmin(7) then V(7)=Xmin(7);
            elseif V(7)>Xmax(7) then V(7)=Xmax(7);
            end
            if V(8)<=Xmin(8) then V(8)=Xmin(8);
            elseif V(8)>Xmax(8) then V(8)=Xmax(8);
            end
            if V(9)<=Xmin(9) then V(9)=Xmin(9);
            elseif V(9)>Xmax(9) then V(9)=Xmax(9);
            end
            if V(10)<=Xmin(10) then V(10)=Xmin(10);
            elseif V(10)>Xmax(10) then V(10)=Xmax(10);
            end
            if V(11)<=Xmin(11) then V(11)=Xmin(11);
            elseif V(11)>Xmax(11) then V(11)=Xmax(11);
            end
            if V(12)<=Xmin(12) then V(12)=Xmin(12);
            elseif V(12)>Xmax(12) then V(12)=Xmax(12);
            end
            if V(13)<=Xmin(13) then V(13)=Xmin(13);
            elseif V(13)>Xmax(13) then V(13)=Xmax(13);
            end
            if V(14)<=Xmin(14) then V(14)=Xmin(14);
            elseif V(14)>Xmax(14) then V(14)=Xmax(14);
            end
            if V(15)<=Xmin(15) then V(15)=Xmin(15);
            elseif V(15)>Xmax(15) then V(15)=Xmax(15);
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
    disp(bM);
    disp(val(bestIndex));
    
//    iterVec=1:1:itermax;
//    plot(iterVec,costVec,2)
endfunction
