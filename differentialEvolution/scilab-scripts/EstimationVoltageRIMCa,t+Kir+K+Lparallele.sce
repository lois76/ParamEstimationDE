/////////////////////////////////////////////////////////////
///////////////     Récupération données      ///////////////
/////////////////////////////////////////////////////////////

a = read("/scilab-scripts/Fig1ARIMCurrentClampTrace.txt",-1,12);
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
    Hdot=zeros(5,1);
    Hdot(1)=(1/pa(26))*(-pa(1)*x(2)*x(3)*(x(1)-pa(5)) - pa(2)*xinf(x(1),pa(10),pa(15))*(x(1)-pa(6)) - pa(3)*x(4)*x(5)*(x(1)-pa(6)) - pa(4)*(x(1)-pa(7)) + I)
    Hdot(2)=(xinf(x(1),pa(8),pa(13))-x(2))/pa(18)
    Hdot(3)=(xinf(x(1),pa(9),pa(14))-x(3))/pa(19)
    Hdot(4)=(xinf(x(1),pa(11),pa(16))-x(4))/pa(20)
    Hdot(5)=(xinf(x(1),pa(12),pa(17))-x(5))/pa(21)
endfunction

//Fonction coût 
function y=fct11(pa)
    c=0;
    condini = [-38; pa(22); pa(23); pa(24); pa(25)]
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


////////////////////////////////////////////////////////
/////////    Estimation de la capacitance C    /////////
////////////////////////////////////////////////////////

function [bM, valBest]=simulation(NP,itermax,F,CR)
    
    D=26;
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.1 0.1 0.1 0.1 20  -100 -90 -90 -90 -90 -90 -90 1  -30 -30 1  -30 0.0001 0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001 0.001];
    Xmax=[50  50  50  50  150 -2   30  -2  -2  -2  -2  -2   30 -1  -1  30 -1  15     15     15     15     0.999 0.999 0.999 0.999 10];
    
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
            if V(17)<=Xmin(17) then V(17)=Xmin(17);
            elseif V(17)>Xmax(17) then V(17)=Xmax(17);
            end
            if V(18)<=Xmin(18) then V(18)=Xmin(18);
            elseif V(18)>Xmax(18) then V(18)=Xmax(18);
            end
            if V(19)<=Xmin(19) then V(19)=Xmin(19);
            elseif V(19)>Xmax(19) then V(19)=Xmax(19);
            end
            if V(20)<=Xmin(20) then V(20)=Xmin(20);
            elseif V(20)>Xmax(20) then V(20)=Xmax(20);
            end
            if V(21)<=Xmin(21) then V(21)=Xmin(21);
            elseif V(21)>Xmax(21) then V(21)=Xmax(21);
            end
            if V(22)<=Xmin(22) then V(22)=Xmin(22);
            elseif V(22)>Xmax(22) then V(22)=Xmax(22);
            end
            if V(23)<=Xmin(23) then V(23)=Xmin(23);
            elseif V(23)>Xmax(23) then V(23)=Xmax(23);
            end
            if V(24)<=Xmin(24) then V(24)=Xmin(24);
            elseif V(24)>Xmax(24) then V(24)=Xmax(24);
            end
            if V(25)<=Xmin(25) then V(25)=Xmin(25);
            elseif V(25)>Xmax(25) then V(25)=Xmax(25);
            end
            if V(26)<=Xmin(26) then V(26)=Xmin(26);
            elseif V(26)>Xmax(26) then V(26)=Xmax(26);
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
    
endfunction





