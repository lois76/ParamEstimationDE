vecV=[-100:10:50]
Inf=[-12.2 -9.13 -6.57 -4.91 -3.57 -2.13 -0.807 0.229 1.46 4.27 7.46 11.8 17.2 21.6 27.1 32.5]

function y=xinf(V,V12,k)
    y=1 ./(1+exp((V12-V) ./k));
endfunction

// Fonction coût pour un modèle de type 2-1
function y=W21(pa)
    e=0;
    for i=1:length(vecV)
        e=e+(Inf(i)-(pa(1)*xinf(vecV(i),pa(7),pa(10))*(vecV(i)-pa(4))+pa(2)*xinf(vecV(i),pa(8),pa(11))*xinf(vecV(i),pa(9),pa(12))*(vecV(i)-pa(5))+pa(3)*(vecV(i)-pa(6))))^2
    end
    y=e
endfunction


/////////////////////////////////////////////////////////////////////////////////////////
//// Estimation de [gCa gK gL ECa EK EL V1/2x1 V1/2x2 V1/2x3 V12x4 kx1 kx2 kx3 kx4] ////
/////////////////////////////////////////////////////////////////////////////////////////

function [valBest]=simulation(NP,itermax,F,CR)
    
    D=12; 
    costVec=zeros(1,itermax);
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0 0 0 20 -90 -80 -90 -90 -90 0 0 -30];
    Xmax=[30 30 30 140 0 30 0 0 0 30 30 0];
    
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
        val(j)=W21(pop(:,j))
    end
    
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
            if V(1)<=0.1 then V(1)=0.1;
            elseif V(1)>30 then V(1)=30;
            end
            if V(2)<=0.1 then V(2)=0.1;
            elseif V(2)>30 then V(2)=30;
            end
            if V(3)<=0.1 then V(3)=0.1;
            elseif V(3)>30 then V(3)=30;
            end
            if V(4)<=20 then V(4)=20;
            elseif V(4)>140 then V(4)=140;
            end
            if V(5)<=-90 then V(5)=-90;
            elseif V(5)>-0.1 then V(5)=-0.1;
            end
            if V(6)<=-80 then V(6)=-80;
            elseif V(6)>30 then V(6)=30;
            end
            if V(7)<=-90 then V(7)=-90;
            elseif V(7)>-2 then V(7)=-2;
            end
            if V(8)<=-90 then V(8)=-90;
            elseif V(8)>-2 then V(8)=-2;
            end
            if V(9)<=-90 then V(9)=-90;
            elseif V(9)>-2 then V(9)=-2;
            end
            if V(10)<=0.1 then V(10)=0.1;
            elseif V(10)>30 then V(10)=30;
            end
            if V(11)<=0.1 then V(11)=0.1;
            elseif V(11)>30 then V(11)=30;
            end
            if V(12)<=-30 then V(12)=-30;
            elseif V(12)>-0.1 then V(12)=-0.1;
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
            tempval=W21(U(:,j));

            if tempval<=val(j) then
                pop(:,j) = U(:,j);
                val(j) = tempval;
            end
        end
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
    
    // Sauvegarde du meilleur individu
    bM = [];
    bM = pop(:,bestIndex);

    valBest=val(bestIndex);
endfunction