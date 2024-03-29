////////////////////////////////////////////////
//// Définition du benchmark HH de type 1-1 ////
////////////////////////////////////////////////

bc=[0.3 
    0.12  
    0.86  
    0.44  
    35.3  
  -77.1  
  -64.23  
  -29.75  
  -74.6   
  -48.9  
  -79.6  
    29.1  
  -15.1  
    28.71  
  -21.0
   0.23
   1.1
   3.7
   0.52
   0.36
   0.02
   0.032]

gCa=bc(1); gKir=bc(2); gK=bc(3); gL=bc(4);
ECa=bc(5); EK=bc(6); EL=bc(7);
V12mCa=bc(8); V12hKir=bc(9); V12mK=bc(10); V12hK=bc(11);
kmCa=bc(12); kKir=bc(13); kmK=bc(14); khK=bc(15);

tmCa=bc(16); tmK=bc(17); thK=bc(18);
mCa0=bc(19); mK0=bc(20); hK0=bc(21);
C=bc(22);

function y=xinf(V,V12,k)
    y=1/(1+exp((V12-V)/k));
endfunction

function [Hdot]=HHbench(t,x,pa)
    Hdot=zeros(4,1);
    Hdot(1)=(1/C)*(-gCa*x(2)*(x(1)-ECa) - gKir*xinf(x(1),V12hKir,kKir)*(x(1)-EK) - gK*x(3)*x(4)*(x(1)-EK) - gL*(x(1)-EL) + I)
    Hdot(2)=(xinf(x(1),V12mCa,kmCa)-x(2))/tmCa
    Hdot(3)=(xinf(x(1),V12mK,kmK)-x(3))/tmK
    Hdot(4)=(xinf(x(1),V12hK,khK)-x(4))/thK
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
    x=ode([-36; mCa0; mK0; hK0],t0,t,HHbench);
    x1=x(1,:);
    plot2d(t,x1,3)
    x1=x1';
    a(:,i+1)=x1;
end


/////////////////////////////////////////////////////////////////////////////////
///////////////    Paramètres déterminés à partir de HillVallEA    //////////////
/////////////////////////////////////////////////////////////////////////////////

par=[0.26296128553453496091 0.11595616381322577937 20.75291714112690755201 0.46226542403717535423 25.02871154943775522383 -86.89227180270813732932 -55.28830786920131146189 -11.44616404922916608200 -59.56141476268292223040 -27.09811870977855363662 -89.99985887251865790404 27.18698184549467100624 -14.73495643454868009314 9.90767134294620355206 -7.97925059809124626042]

//par=[0.2925687
//    0.1144681  
//    0.2561783  
//    0.4477050  
//    41.453044  
//  -70.557332  
//  -66.482265  
//  -28.3064    
//  -50.256273  
//  -84.004184  
//  -79.199446  
//    28.956379  
//  -21.358954  
//    29.017845  
//  -16.099689 
//]

/////////////////////////////////////////////////////////////
///////////////    Fonction coût algorithme    //////////////
/////////////////////////////////////////////////////////////


function [Hdot]=HH(t,x,pa)
    Hdot=zeros(4,1);
    Hdot(1)=(1/pa(7))*(-par(1)*x(2)*(x(1)-par(5)) - par(2)*xinf(x(1),par(9),par(13))*(x(1)-par(6)) - par(3)*x(3)*x(4)*(x(1)-par(6)) - par(4)*(x(1)-par(7)) + I)
    Hdot(2)=(xinf(x(1),par(8),par(12))-x(2))/pa(1)
    Hdot(3)=(xinf(x(1),par(10),par(14))-x(3))/pa(2)
    Hdot(4)=(xinf(x(1),par(11),par(15))-x(4))/pa(3)
endfunction

//Fonction coût 

function y=fct(pa)
    c=0;
    condini = [-36; pa(4); pa(5); pa(6)]
    for i=1:11
        I=stim(i);
        x=ode(condini,t0,t,HH); 
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



