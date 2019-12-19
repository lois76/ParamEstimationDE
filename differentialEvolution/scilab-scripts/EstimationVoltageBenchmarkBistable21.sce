/////////////////////////////////////////////////////////////
///////////////     Récupération données      ///////////////
/////////////////////////////////////////////////////////////

gCa=0.52; gK=5; gL=0.8;
ECa=70; EK=-88; EL=-73;
V12x1=-35; V12x2=-26; V12x3=-62;
kx1=5; kx2=14; kx3=-10;
tx1=1.6; tx2=2.7; tx3=3.8; 
x10=0.9; x20=0.61; x30=0.21;
C=0.16;

function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

function [Hdot]=HH21ref(t,x)
    Hdot=zeros(4,1);  //définition de la taille
    Hdot(1)=(1/C)*(-gCa*x(2)*(x(1)-ECa)-gK*x(3)*x(4)*(x(1)-EK)-gL*(x(1)-EL)+I)
    Hdot(2)=(xinf(x(1),V12x1,kx1)-x(2))/tx1
    Hdot(3)=(xinf(x(1),V12x2,kx2)-x(3))/tx2
    Hdot(4)=(xinf(x(1),V12x3,kx3)-x(4))/tx3
endfunction

t=linspace(0,40,20000);
t0=0;
a=zeros(length(t),11);
stim=[-15:5:35];

for i=1:11
    I=stim(i);
    x=ode([-46; x10; x20; x30],t0,t,HH21ref);
    x1=x(1,:);
    x1=x1';
    a(:,i)=x1;
//    plot2d(t,x1,2)
end





/////////////////////////////////////////////////////////////
///////////////    Fonction coût algorithme    //////////////
/////////////////////////////////////////////////////////////

function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

function [Hdot]=HH21est(t,x,pa)
    Hdot=zeros(4,1);
    Hdot(1)=(1/pa(19))*(-pa(1)*x(2)*(x(1)-pa(4))-pa(2)*x(3)*x(4)*(x(1)-pa(5))-pa(3)*(x(1)-pa(6))+I)
    Hdot(2)=(xinf(x(1),pa(7),pa(10))-x(2))/pa(13)
    Hdot(3)=(xinf(x(1),pa(8),pa(11))-x(3))/pa(14)
    Hdot(4)=(xinf(x(1),pa(9),pa(12))-x(4))/pa(15)
endfunction

//Fonction coût 

t0=0;
function y=fct21(pa)
    c=0;
    condini = [-46; pa(16); pa(17); pa(18)]
    for i=1:11
        I=stim(i);
        x=ode(condini,t0,t,HH21est); 
        Vest=x(1,:);
        for k=1:length(t)
            c=c+(Vest(k)-a(k,i))*(Vest(k)-a(k,i))
        end
    end
    y=c/length(t);
endfunction



////////////////////////////////////////////////////////
/////////    Estimation de la capacitance C    /////////
////////////////////////////////////////////////////////

function [bM]=SS(NP,itermax,F,CR)
    
    D=19;
//    costVec=zeros(1,itermax);
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.1 0.1 0.1 20  -90 -80 -90 -90 -90  1  1  -30 0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001];
    Xmax=[50  50  50  150 -2  30  -2  -2  -2   30 30 -1  15     15     15     0.999 0.999 0.999 10];
    
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
        val(j)=fct21(pop(:,j))
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
            tempval=fct21(U(:,j));

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
//    disp(bestIndex);
    
    // Sauvegarde du meilleur individu
    bM = [];
    bM = pop(:,bestIndex);
    
    disp(val);
    disp(bM);
    disp(val(bestIndex));
    
//    iterVec=1:1:itermax;
//    plot(iterVec,costVec,2)
endfunction

[bM]=SS(180,1100,0.5,0.85)




