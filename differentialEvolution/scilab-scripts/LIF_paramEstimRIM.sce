/////////////////////////////////////////////////////
///////////////     Voltage Data      ///////////////
/////////////////////////////////////////////////////

//a = read("/home/naudin/Documents/FichierScilab/EstimationRIM/Fig1ARIMCurrentClampTrace.txt",-1,12);
a = read("/scilab-scripts/Fig1ARIMCurrentClampTrace.txt",-1,12);
A=a(2489:14988,2:$)*1000;
t=linspace(0,50,12500);
t0=0;
stim=[-15:5:35];

/////////////////////////////////////////////////////////////
///////////////    Fonction coût algorithme    //////////////
/////////////////////////////////////////////////////////////

function [Hdot]=HH(t,x,pa)
    Hdot(1)=(pa(1)/pa(2))*(x(1)/pa(1) + I)
endfunction

%ODEOPTIONS=[1,0,0,%inf,0,2,20000,12,5,0,-1,-1];

function y=fct11(pa)
    c=0;
    condini = -38
    for i=1:11
        I=stim(i);
        x=ode(condini,t0,t,HH); 
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-A(k,i))*(V(k)-A(k,i))
        end
    end
    y=c/length(t);
endfunction

function [bM, valBest]=simulation(NP,itermax,F,CR)
    
    D=2;//nombre de paramètres à estimer
    pop=zeros(D,NP);//population : matrice de taille D*NP où NP est le nombre d'individus de ta population

    ////////////////////////////////////////////////////////////////////////
    //// Vecteurs de contraintes des paramètres : borne minimum/maximum ////
    ////////////////////////////////////////////////////////////////////////

    Xmin=[1   1];
    Xmax=[100 100];
    
    ///////////////////////////////////////////////////
    //// Initialisation aléatoire de ma population ////
    ///////////////////////////////////////////////////

    for j=1:NP
        for i=1:D
            pop(i,j)=Xmin(i)+(Xmax(i)-Xmin(i))*rand();
        end
    end
    
    //////////////////////////////////////////////////////////////
    //// Évaluation du meilleur individu après initialisation ////
    //////////////////////////////////////////////////////////////
    
    val=zeros(NP,1); // tableau avec le coût de chacun des individus, initialisé à 0
    
    // Evaluation de la fonction coût pour chacun des individus de la population
    for j=1:NP
        val(j)=fct11(pop(:,j))
    end
    disp(val)
    
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
            
            // ======== Contraintes ========
            for i=1:D
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
            tempval=fct11(U(:,j));

            if tempval<=val(j) then
                pop(:,j) = U(:,j);
                val(j) = tempval;
            end
        end
        disp(iter) // c'est juste pour voir l'avancée de l'algorithme, à quelle itération on est...
        iter = iter + 1;
    end  //fin de la boucle while
    
    // Détermination de l'indice du meilleur individu
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
    valBest=val(bestIndex);
    
    // Sauvegarde du meilleur individu
    bM = [];
    bM = pop(:,bestIndex);
    
    disp(val);
    disp(pop)
endfunction

//simulation(20,10,0.5,0.9)
