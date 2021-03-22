/////////////////////////////////////////////////////////////
///////////////     Récupération données      ///////////////
/////////////////////////////////////////////////////////////

//Voltage
a = read("/scilab-scripts/Fig1ARIMCurrentClampTrace.txt",-1,12);
//a = read("/home/naudin/Documents/FichierScilab/Fourre tout/Fig1A_RIMCurrentClampTrace.txt",-1,12);
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

function [xdot]=HR(t,x,pa)
//    xdot=zeros(3,1);
    xdot(1) = x(2) - pa(1)*x(1)*x(1)*x(1) + pa(2)*x(1)*x(1) - x(3) + I
    xdot(2) = pa(3) - pa(4)*x(1)*x(1) - x(2)
    xdot(3) = pa(5)*(pa(6)*(x(1)-x0)-x(3))
endfunction

x0=-38

//pa(1)=1; pa(2)=3; pa(3)=1; pa(4)=5.0; pa(5)=0.006; pa(6)=4.0; pa(7)=-1.6;
//I=2
//x=ode([-1.2; -12; 0],t0,t,HRbis);
//V=x(1,:);
//plot2d(t,V,2)

%ODEOPTIONS=[1,0,0,%inf,0,2,20000,12,5,0,-1,-1];

//Fonction coût 
function y=fct11(pa)
    c=0;
    condini = [-38; pa(7); pa(8)]
    for i=1:11
        I=stim(i);
        x=ode(condini,t0,t,HR); 
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-A(k,i))*(V(k)-A(k,i))
        end
    end
    y=c/length(t);
endfunction

/////////////////////////////////////////////
/////////    Parameter estimation   /////////
/////////////////////////////////////////////

function [bM, valBest]=simulation(NP,itermax,F,CR)
    
    D=8;
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.001 0.001 0.001 0.001 0.0001 0.001   -50 -50];
    Xmax=[50    50    50    50    1      50  50  50];
    
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
        val(j)=fct11(pop(:,j))
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
    disp(val)
endfunction



