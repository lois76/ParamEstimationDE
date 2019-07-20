/////////////////////////////////////////////////////////////
///////////////     Simulation Benchmark      ///////////////
/////////////////////////////////////////////////////////////

function y=xinf(V,V12,k)
    y=1 ./(1+exp((V12-V) ./k));
endfunction

//Paramètres du benchmark
gCa=0.4;  gK=0.9; gL=0.24;
ECa=126.5; EK=-86; EL=-49;
V12x1=-4.3;  V12x2=-13.01;  V12x3=-51.7; V12x4=-82.2;
kx1=22.8; kx2=-29.9; kx3=10.8; kx4=-3.3;

//
function y=Iinf(VH)
    y=gCa.*xinf(VH,V12x1,kx1).*xinf(VH,V12x2,kx2).*(VH-ECa)+gK.*xinf(VH,V12x3,kx3).*xinf(VH,V12x4,kx4).*(VH-EK)+gL.*(VH-EL)
endfunction

vecV=[-100:10:50];


//pa=[25.31913545036314 7.089103872041187 0.4953453153730928 119.20016463657245 -47.61527118817742 -17.625711243895452 -6.455069029894645 -75.43521676361777 -24.470271226478673 -85.23875139930527 9.751113286530472 -17.22756685182511 2.5749795352730755 -23.67912116804757]

//Fonction coût 
function y=W(pa)
    e=0;
    for i=1:length(vecV)
        e=e+(Iinf(vecV(i))-(pa(1)*xinf(vecV(i),pa(7),pa(11))*xinf(vecV(i),pa(8),pa(12))*(vecV(i)-pa(4))+pa(2)*xinf(vecV(i),pa(9),pa(13))*xinf(vecV(i),pa(10),pa(14))*(vecV(i)-pa(5))+pa(3)*(vecV(i)-pa(6))))^2
    end
    y=e/length(vecV)
endfunction

//disp(W(pa))


////////////////////////////////////////////////////////////////////////////////////////
//// Estimation de [gCa gK gL ECa EK EL V1/2x1 V1/2x2 V1/2x3 V12x4 kx1 kx2 kx3 kx4] ////
////////////////////////////////////////////////////////////////////////////////////////

function [bM, valBest]=simulation(NP,itermax,F,CR)
    
    D=14; 
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.1 0.1 0.1 20  -90 -80 -90 -90 -90 -90 1  -30 1  -30];
    Xmax=[30  30  30  150 -2  30  -2  -2  -2  -2  30 -1  30 -1];
    
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
    
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
    
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
        iter = iter + 1;
        bestIndex=1;
        for b=2:NP
            if val(b)<val(bestIndex) then bestIndex=b; end
        end
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

    valBest = val(bestIndex);
    
//    iterVec=1:1:itermax;
//    plot(iterVec,costVec,2)
endfunction
