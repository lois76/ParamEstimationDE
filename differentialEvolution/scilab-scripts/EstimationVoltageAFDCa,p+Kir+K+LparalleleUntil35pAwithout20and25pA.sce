/////////////////////////////////////////////////////////////
///////////////     Récupération données      ///////////////
/////////////////////////////////////////////////////////////

a = read("/scilab-scripts/AFDnewDataSecondRecordingsNumberPointsDividedBy4.txt",-1,11);
//a = read("/home/naudin/Documents/article-2/AFD under Extreme Stimulation/Second Recordings/AFDnewDataSecondRecordingsNumberPointsDividedBy4.txt",-1,11);
t=linspace(0,50,12501);
//A=a(:,1:$);
t0=0;
stim1=[-15:5:15];
A1=a(:,1:length(stim1));
stim2=[30 35]
A2=a(:,10:11);
//stim=[stim1 stim2]

//for i=[1:1:size(A2,'c')]
//    plot2d(t,A2(:,i),3)
//end

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

//Function that computes the standard deviation
function y=sigma(v)
    s=0;
    moy=mean(v);
    for i=1:length(v)
        s=s+(v(i)-moy)^2
    end
    y=sqrt(s/(length(v)-1));
endfunction

//Noise level (standard deviation) for each I
dev1=[]
for i=1:length(stim1)
    dev1=[dev1 sigma(A1(10000:$,i))]
end
dev2=[sigma(A2(10000:$,1)) sigma(A2(10000:$,2))]

%ODEOPTIONS=[1,0,0,%inf,0,2,20000,12,5,0,-1,-1];
//Cost function voltage
function y=fct11(pa)
    tmp=0;
    condini = [-76; pa(19); pa(20); pa(21)]
    for i=1:length(stim1)
        c=0;
        I=stim1(i);
        x=ode(condini,t0,t,HH11); 
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-A1(k,i))*(V(k)-A1(k,i));
        end
        c=sqrt(c/length(t))/dev1(i);
        tmp=tmp+c;
    end
    for i=1:length(stim2)
        c=0;
        I=stim2(i)
        x=ode(condini,t0,t,HH11); 
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-A2(k,i))*(V(k)-A2(k,i));
        end
        c=sqrt(c/length(t))/dev2(i);
        tmp=tmp+c;
    end
    y=tmp/(length(stim1)+length(stim2));
endfunction

////////////////////////////////////////////////////////
/////////    Estimation de la capacitance C    /////////
////////////////////////////////////////////////////////

function [bM, valBest]=simulation(NP,itermax,F,CR)
    
    D=22;//nombre de paramètres du modèle
    pop=zeros(D,NP);//matrice de taille DxNP (i.e. D lignes et NP colonnes, où chaque colonne contient un vecteur de D paramètres) 

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.0001 0.0001 0.0001 0.0001 20  -100 -90 -90 -90 -90 -90 1  -30 1  -30 0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001];
    Xmax=[50     50     50     50     150 -2   30  -2  -2  -2  -2  30 -1  30 -1  20     20     20     0.999 0.999 0.999 10];
    
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
    
    //Évaluation de chacun des individus
    for j=1:NP
        val(j)=fct11(pop(:,j))
    end
    
    //On détermine le meilleur individu de la population, i.e. celui qui a la plus petite fonction coût
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
            
            // ======== Contraintes sur les paramètres ========
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
        disp(iter);//juste une commande pour voir où tu en es dans les itérations
        iter = iter + 1;

//	//détermination du meilleur individu (pourquoi ? c'est une bonne question -_-')
//        bestIndex=1;
//        for b=2:NP
//            if val(b)<val(bestIndex) then bestIndex=b; end
//        end
    end  //fin de la boucle while
    


    // Détermination de l'indice du meilleur individu, noté bestIndex
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
    valBest=val(bestIndex);//c'est la valeur coût du meilleur individu (i.e. celui avec la plus petite fonction coût) 
    
    // Sauvegarde du meilleur individu
    bM = [];
    bM = pop(:,bestIndex);
    
endfunction

//[bM, valBest]=simulation(10,5,0.5,0.9)



