////////////////////////////////////////////////
//// Définition du benchmark HH de type 1-1 ////
////////////////////////////////////////////////

bc=[0.2 
    0.4 
    0.3
    104
    -18.4
    -63.1
    -23.8
    -37.2
    -25.1
    18.4
    -27.8
    16.3
    0.6
    5
    2.3
    0.3
    0.2
    0.1
    0.04]

gCa=bc(1);  gK=bc(2); gL=bc(3);
ECa=bc(4); EK=bc(5); EL=bc(6);
V12mCa=bc(7); V12hCa=bc(8); V12mK=bc(9);
kmCa=bc(10); khCa=bc(11); kmK=bc(12);
tmCa=bc(13); thCa=bc(14); tmK=bc(15);
mCa0=bc(16); hCa0=bc(17); mK0=bc(18);
C=bc(19);

function y=xinf(V,V12,k)
    y=1/(1+exp((V12-V)/k));
endfunction

function [Hdot]=HHbench(t,x,pa)
    Hdot=zeros(4,1);
    Hdot(1)=(1/C)*(-gCa*x(2)*x(3)*(x(1)-ECa) - gK*x(4)*(x(1)-EK) - gL*(x(1)-EL) + I)
    Hdot(2)=(xinf(x(1),V12mCa,kmCa)-x(2))/tmCa
    Hdot(3)=(xinf(x(1),V12hCa,khCa)-x(3))/thCa
    Hdot(4)=(xinf(x(1),V12mK,kmK)-x(4))/tmK
endfunction

/// Stimulations appliquées au neurone ///
stim=[-15:5:35];

/// Construction des solutions ///
t0=0;
t=linspace(0,50,12500);
t=t';
a=zeros(length(t),12);
a(:,1)=t;
for i=1:11
    I=stim(i);
    x=ode([-40;mCa0;hCa0;mK0],t0,t,HHbench);
    x1=x(1,:);
//    plot2d(t,x1,3)
    x1=x1';
    a(:,i+1)=x1;
end

//write('/home/loisse/Documents/FichierScilab/Article 2/CurrentClampBenchmarkLinear1Ca,t+K,p+L.txt', a)

/////////////////////////////////////////////////////////////////////////////////
///////////////    Paramètres déterminés à partir de HillVallEA    //////////////
/////////////////////////////////////////////////////////////////////////////////

//par=[2.38664969875664789711 0.43734299253905523086 0.28283396456824699827 76.20810707090238622641 -18.99141337004144247658 -58.67142768132298868977 -24.96397259325666695418 -74.51408155882776895851 -11.96464196615497144194 10.25971902778186795047 -14.20947412306290047468 13.59353177812722002216]

par=[0.3621152  
    0.4047906  
    0.2997707  
    49.805558  
  -16.82725   
  -63.025241  
  -23.66706   
  -44.294987  
  -25.052879  
    16.428948  
  -30.834623  
    17.402578]



/////////////////////////////////////////////////////////////
///////////////    Fonction coût algorithme    //////////////
/////////////////////////////////////////////////////////////

function [Hdot]=HH12(t,x,pa)
    Hdot=zeros(4,1);
    Hdot(1)=(1/pa(7))*(-par(1)*x(2)*x(3)*(x(1)-par(4)) - par(2)*x(4)*(x(1)-par(5)) - par(3)*(x(1)-par(6)) + I)
    Hdot(2)=(xinf(x(1),par(7),par(10))-x(2))/pa(1)
    Hdot(3)=(xinf(x(1),par(8),par(11))-x(3))/pa(2)
    Hdot(4)=(xinf(x(1),par(9),par(12))-x(4))/pa(3)
endfunction

//Fonction coût 

function y=fct(pa)
    c=0;
    condini = [-40; pa(4); pa(5); pa(6)]
    for i=1:11
        I=stim(i);
        x=ode(condini,t0,t,HH12); 
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-a(k,i+1))*(V(k)-a(k,i+1))
        end
    end
    y=c/length(t);
endfunction

//////////////////////////////////////
////     Parameter estimation     ////
//////////////////////////////////////

function [bM, valBest]=simulation(NP,itermax,F,CR)

    D=7; 
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////
    
    Xmin=[0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001];
    Xmax=[15     15     15     0.999 0.999 0.999 10];

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
        val(j)=fct(pop(:,j))
    end
    
    disp(val)
    
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
            r1=j; r2=j; r3=j;
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
            tempval=fct(U(:,j));

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
endfunction



