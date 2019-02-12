currentDate = getdate()
rand('seed',getdate('s')+currentDate(10))

//////////////////////////////////////////////////////////
///////////////     Simulation modèle      ///////////////
//////////////////////////////////////////////////////////

function y=xinf(V,V12,k)
    y=1 ./(1+exp((V12-V) ./k));
endfunction

ENa=120; EK=-12; EL=13;
gNa=2; gK=5; gL=0.8;

function y=Iinf(VH)
    y=gNa.*xinf(VH,-35,5).*xinf(VH,-52,-7).*(VH-ENa)+gK.*xinf(VH,-26,14).*xinf(VH,-62,-10).*(VH-EK)+gL.*(VH-EL)
endfunction

// Définition des données à partir desquelles on va travailler
vecV=linspace(-100,50,16)

Inf=[];
for j=1:length(vecV)
    Inf(j)=Iinf(vecV(j));
end

function y=W(pa)
    e=0;
    for i=1:length(vecV)
        e=e+(Inf(i)-(pa(1)*xinf(vecV(i),pa(7),pa(11))*xinf(vecV(i),pa(8),pa(12))*(vecV(i)-pa(4))+pa(2)*xinf(vecV(i),pa(9),pa(13))*xinf(vecV(i),pa(10),pa(14))*(vecV(i)-pa(5))+pa(3)*(vecV(i)-pa(6))))^2
    end
    y=e
endfunction


/////////////////////////////////////////////////////////////////////////////////////////
//// Estimation de [gCa gK gL ECa EK EL V1/2x1 V1/2x2 V1/2x3 V12x4 kx1 kx2 kx3 kx4] ////
/////////////////////////////////////////////////////////////////////////////////////////

function [valBest]=SS(NP,itermax,F,CR)
    D=14; 
    costVec=zeros(1,itermax);
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0 0 0 80 -50 0 -70 -70 -70 -70 0 -15 0 -15];
    Xmax=[10 10 10 140 0 20 0 0 0 0 15 0 15 0];
    
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
            elseif V(1)>10 then V(1)=10;
            end
            if V(2)<=0.1 then V(2)=0.1;
            elseif V(2)>10 then V(2)=10;
            end
            if V(3)<=0.1 then V(3)=0.1;
            elseif V(3)>10 then V(3)=10;
            end
            if V(4)<=80 then V(4)=80;
            elseif V(4)>140 then V(4)=140;
            end
            if V(5)<=-50 then V(5)=-50;
            elseif V(5)>-0.1 then V(5)=-0.1;
            end
            if V(6)<=0.1 then V(6)=0.1;
            elseif V(6)>20 then V(6)=20;
            end
            if V(7)<=-70 then V(7)=-70;
            elseif V(7)>-2 then V(7)=-2;
            end
            if V(8)<=-70 then V(8)=-70;
            elseif V(8)>-2 then V(8)=-2;
            end
            if V(9)<=-70 then V(9)=-70;
            elseif V(9)>-2 then V(9)=-2;
            end
            if V(10)<=-70 then V(10)=-70;
            elseif V(10)>-2 then V(10)=-2;
            end
            if V(11)<=0.1 then V(11)=0.1;
            elseif V(11)>15 then V(11)=15;
            end
            if V(12)<=-15 then V(12)=-15;
            elseif V(12)>-0.1 then V(12)=-0.1;
            end
            if V(13)<=0.1 then V(13)=0.1;
            elseif V(13)>15 then V(13)=15;
            end
            if V(14)<=-15 then V(14)=-15;
            elseif V(14)>-0.1 then V(14)=-0.1;
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
//        disp(iter)
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
    
//    disp(val);
//    disp(bM);
    disp(val(bestIndex));
    
    
//    iterVec=1:1:itermax;
//    plot(iterVec,costVec,2)
endfunction

args=sciargs();
[valBest]=SS(strtod(args(6)),strtod(args(7)),strtod(args(8)),strtod(args(9)));

csvWrite(valBest,"/Result/result.csv")
