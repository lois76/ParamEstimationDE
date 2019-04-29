// Récupération des données expérimentales sous forme de matrice
a = read("/scilab-scripts/DataRIMVoltageClamp.txt",-1,17);
A=a(239:1488,2:$)*1E12;
t=linspace(0,50,1250);
//write('/home/loisse/Documents/FichierScilab/EstimationSteadyStateCurent/EstimationSSLinear/A.txt',a);

////////////////////////////////////////////////////////////////
//// Estimation de [tx1 tx2 tx3 tx4 x1(0) x2(0) x3(0) x4(0) ////
////////////////////////////////////////////////////////////////

////Premier jeu de paramètres : fonction coût = 0.3658
//gCa=0.175; gK=29.65; gL=0.357;
//ECa=103.15; EK=-5.996; EL=-66.719;
//V12x1=-26.50; V12x2=-17.29; V12x3=-66.21;
//kx1=11.60; kx2=8.89; kx3=-9.54;

////Deuxième jeu de paramètres : fonction coût = 0.3650
gCa=0.175; gK=30; gL=0.366;
ECa=113.32; EK=-38.646; EL=-68.808;
V12x1=-34.411; V12x2=-37.415; V12x3=-69.735;
kx1=16.995; kx2=5.225; kx3=-5.658;

////Trosième jeu de paramètres : fonction coût = 0.3686
//gCa=0.185; gK=0.255; gL=0.369;
//ECa=113.964; EK=-33.134; EL=-68.44;
//V12x1=-32.137; V12x2=-61.4288; V12x3=-41.539;
//kx1=14.646; kx2=8.5508; kx3=-24.77;

function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

function y=x(VH,V12,k,tx,x0,t)
    y=xinf(VH,V12,k)+(x0-xinf(VH,V12,k))*exp(-t/tx)
endfunction

function y=Iest(VH,pa,t)
    y=gCa.*x(VH,V12x1,kx1,pa(1),pa(4),t).*(VH-ECa) + gK.*x(VH,V12x2,kx2,pa(2),pa(5),t).*x(VH,V12x3,kx3,pa(3),pa(6),t).*(VH-EK) + gL.*(VH-EL)
endfunction

VH=[-100:10:50];

function y=W(pa)
    c=0;
    for i=1:length(VH)
        for k=1:length(t)
            c=c+(Iest(VH(i),pa,t(k))-A(k,i)).^2
        end
    end
    y=c/length(t);
endfunction



////////////////////////////////////////////////////////
//// Définition Algorithme Évolution Différentielle ////
////////////////////////////////////////////////////////

function [bestMember, valBest]=simulation(NP,itermax,F,CR)
    
    D=6; 
//    costVec=zeros(1,itermax);
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.01 0.01 0.01 0.001 0.001 0.001];
    Xmax=[15 15 15 0.999 0.999 0.999];

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

    val = zeros(NP,1); // tableau avec le coût de chacun des individus

    // --------- Evalutation du premier individu après initialisation --------

    for j=1:NP
        val(j)=W(pop(:,j));
    end 
    
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
    costVec(1)=val(bestIndex);

    //// =============== ÉTAPE SUIVANTE ================
    iter=1; // nombre d'itération
    B=zeros(D,NP); // Vecteur intermédiaire perturbé (mutation + crossover)
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
            // Condition pour que tous les termes soient différents de 0
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
            // ======== Crossover ========
            for i=1:D
                if rand()<CR then
                    B(i,j)=V(i);
                else
                    B(i,j)=pop(i,j);
                end
            end
        end 
        // ====== Fin de construction de la matrice U =====

        // ======== Sélection ========
        for j=1:NP
            tempval=W(B(:,j));

            if tempval<=val(j) then
                pop(:,j) = B(:,j);
                val(j) = tempval;
            end
        end
        iter = iter + 1;
        bestIndex=1;
        for b=2:NP
            if val(b)<val(bestIndex) then bestIndex=b; end
        end
        costVec(iter)=val(bestIndex);
    end //fin de la boucle while

    
    // Détermination de l'inidice du meilleur individu
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end


    // Sauvegarde du meilleur individu
    bestMember = [];
    bestMember = pop(:,bestIndex);
    valBest = val(bestIndex);
//    iterVec=1:1:itermax
//    plot(iterVec,costVec,2)

endfunction