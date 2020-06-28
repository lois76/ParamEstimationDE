/////////////////////////////////////////////////////////////
///////////////     Récupération données      ///////////////
/////////////////////////////////////////////////////////////

a = read("/scilab-scripts/Fig 1A_AFD Current-Clamp Trace.txt",-1,12);
A=a(2489:14988,2:$)*1000;
t=linspace(0,50,12500);
t0=0;
stim=[-15:5:35];

/////////////////////////////////////////////////////////////
///////////////    Fonction coût algorithme    //////////////
/////////////////////////////////////////////////////////////

function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

function [Hdot]=HH11(t,x,pa)
    Hdot=zeros(4,1);
    Hdot(1)=(1/pa(22))*(-pa(1)*x(2)*(x(1)-pa(5)) - pa(2)*xinf(x(1),pa(9),pa(13))*(x(1)-pa(6)) - pa(3)*x(3)*x(4)*(x(1)-pa(6)) - pa(4)*(x(1)-pa(7)) + I)
    Hdot(2)=(xinf(x(1),pa(8),pa(12))-x(2))/pa(16)
    Hdot(3)=(xinf(x(1),pa(10),pa(14))-x(3))/pa(17)
    Hdot(4)=(xinf(x(1),pa(11),pa(15))-x(4))/pa(18)
endfunction

//Fonction coût 
function y=fct11(pa)
    c=0;
    condini = [-78; pa(19); pa(20); pa(21)]
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


//////////////////////////////////////////
/////////    Estimation Param    /////////
//////////////////////////////////////////

Fl=0.1; Fu=0.9;
t1=0.1; t2=0.1;


function [bM, valBest, costVec, vecF, vecCR]=simulation(NP,itermax)
    
    D=24; // 22+2 params pour F et CR
    costVec=zeros(1,itermax);
    pop=zeros(D,NP);
    
    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////


    Xmin=[0.0001 0.0001 0.0001 0.0001 20  -100 -90 -90 -90 -90 -90 1  -30 1  -30 0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001 0.1 0];
    Xmax=[50     50     50     50     150 -2   30  -2  -2  -2  -2  30 -1  30 -1  15     15     15     0.999 0.999 0.999 10    1   1];
    
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
        val(j)=fct11(pop(:,j))
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
            tempval=fct11(U(:,j));

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
//    disp(vecF);
    
//    iterVec=1:1:itermax;
//    plot(iterVec,costVec,2)
endfunction

