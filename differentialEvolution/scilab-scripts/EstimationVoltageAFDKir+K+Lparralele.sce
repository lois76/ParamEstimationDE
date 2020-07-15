/////////////////////////////////////////////////////////////
///////////////     Récupération données      ///////////////
/////////////////////////////////////////////////////////////

a = read("/home/loisse/Documents/FichierScilab/EstimationAFD/Fig1A_AFDCurrentClampTrace.txt",-1,12);
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

function [Hdot]=HH21(t,x,pa)
    Hdot=zeros(3,1);
    Hdot(1)=(1/pa(16))*(-pa(1)*xinf(x(1),pa(6),pa(9))*(x(1)-pa(4))-pa(2)*x(2)*x(3)*(x(1)-pa(4))-pa(3)*(x(1)-pa(5))+I)
    Hdot(2)=(xinf(x(1),pa(7),pa(10))-x(2))/pa(12)
    Hdot(3)=(xinf(x(1),pa(8),pa(11))-x(3))/pa(13)
endfunction

//Fonction qui calcule l'écart type (standard deviation)
function y=sigma(v)
    s=0;
    moy=mean(v);
    for i=1:length(v)
        s=s+(v(i)-moy)^2
    end
    y=sqrt(s/length(v));
endfunction

dev=[]
for i=1:11
    dev=[dev sigma(A(7000:$,i))]
end

//Fonction coût 
stim=[-15:5:35];
t0=0;
function y=fct11(pa)
    c=0;
    condini = [-78; pa(14); pa(15)]
    for i=1:11
        I=stim(i);
        x=ode(condini,t0,t,HH21); 
        V=x(1,:);
        for k=1:length(t)
            c=c+((V(k)-A(k,i))/dev(i))*((V(k)-A(k,i))/dev(i));
        end
    end
    y=c/length(t);
endfunction


////////////////////////////////////////////////////////
/////////    Estimation de la capacitance C    /////////
////////////////////////////////////////////////////////

function [bM, valBest]=simulation(NP,itermax,F,CR)
    
    D=16;
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.1 0.1 0.1 -100 -90 -100 -90 -90 -30  1  -30  0.0001 0.0001 0.001 0.001 0.001];
    Xmax=[50  50  50  -2   30  -2   -2  -2  -1   30 -1   15     15     0.999 0.999 10];

    /////////////////////////////////////////
    //// Initialisation de ma population ////
    /////////////////////////////////////////

    for j=1:NP
        for i=1:D
            pop(i,j)=Xmin(i)+(Xmax(i)-Xmin(i))*rand();
        end
    end
    
    disp(pop)
    
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
            if V(15)<=Xmin(15) then V(15)=Xmin(15);
            elseif V(15)>Xmax(15) then V(15)=Xmax(15);
            end
            if V(16)<=Xmin(16) then V(16)=Xmin(16);
            elseif V(16)>Xmax(16) then V(16)=Xmax(16);
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
endfunction
