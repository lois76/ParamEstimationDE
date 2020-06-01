/////////////////////////////////////////////////////////////////////////////
///////////////     Récupérations données expérimentales      ///////////////
/////////////////////////////////////////////////////////////////////////////

vecV=[-100:10:50]
Inf=[-12.2 -9.13 -6.57 -4.91 -3.57 -2.13 -0.807 0.229 1.46 4.27 7.46 11.8 17.2 21.6 27.1 32.5]

//plot(vecV,Inf,'o')

/////////////////////////////////////////////////////////////
///////////////    Fonction coût algorithme    //////////////
/////////////////////////////////////////////////////////////

function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

function y=W(pa)
    e=0;
    for i=1:length(vecV)
        e=e+(Inf(i)-(pa(1).*xinf(vecV(i),pa(8),pa(12)).*(vecV(i)-pa(5)) + pa(2).*xinf(vecV(i),pa(9),pa(13)).*(vecV(i)-pa(6)) + pa(3).*xinf(vecV(i),pa(10),pa(14)).*xinf(vecV(i),pa(11),pa(15)).*(vecV(i)-pa(6)) + pa(4).*(vecV(i)-pa(7))))^2
    end
    y=e/length(vecV) 
endfunction


////////////////////////////////////////////////////////
/////////    Estimation de la capacitance C    /////////
////////////////////////////////////////////////////////

Fl=0.1; Fu=0.9;
t1=0.1; t2=0.1;

// Ici, F et CR sont juste les valeurs initiales des 
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
    
    ///////////////////////////////////////////////////////////////////
    //// Évaluation du meilleur individu de ma première génération ////
    ///////////////////////////////////////////////////////////////////
    
    val=zeros(NP,1); // tableau avec le coût de chacun des individus
    
    for j=1:NP
        val(j)=W(pop(:,j))
    end
    
    disp(val);
    
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
//    costVec(1)=val(bestIndex);
    
    ////////////////////////
    //// Étape suivante ////
    ////////////////////////
     
    iter=1; // nombre d'itération
    U=zeros(D,NP); // Vecteur intermédiaire perturbé (mutation + crossover)
    tempval=0;
    
    ////////////////////////////////
    //// Initialisation F et CR ////
    ////////////////////////////////
    
    F=0.5; 
    CR=0.9;
    vecF=zeros(1,itermax);
    vecCR=zeros(1,itermax);
    
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
            
            // ======== Adaptation de F et CR =======
            if rand()<t1 then
                F=Fl+rand()*Fu
            end
            vecF(iter)=F;
            
            if rand()<t2 then
                CR=rand()
            end
            vecCR(iter)=CR;

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
//        costVec(iter)=val(bestIndex);
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


    
