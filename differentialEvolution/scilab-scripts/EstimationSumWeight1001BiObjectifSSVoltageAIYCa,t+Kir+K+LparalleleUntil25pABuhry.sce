//////////////////////////////////////////////////////////
///////////////     Experimental data      ///////////////
//////////////////////////////////////////////////////////

//a = read("/home/naudin/Documents/FichierScilab/Fourre tout/Fig1A_AIYCurrentClampTrace2.txt",-1,12);
a = read("/scilab-scripts/Fig 1A_AIY Current-Clamp Trace.txt",-1,12);
A=a(2489:14988,2:10)*1000;
t=linspace(0,50,12500);
t0=0;
stim=[-15:5:25];

//Steady-state current
vecV=[-120:10:50];
Inf=[-13.1 -10.4 -7.92 -5.89 -4.11 -2.69 -1.02 0.0211 1.17 3.1 7.32 14.2 22.4 31.5 43.2 54.5 69.5 82.4];
InfSD=[2.88 2.55 1.47 1.31 1.04 0.809 0.7 0.658 0.638 0.889 1.94 3.5 5.36 7.63 10.6 13.3 16 17.9]

//////////////////////////////////////////////////
///////////////    Cost function    //////////////
//////////////////////////////////////////////////

//Boltzmann function
function y=xinf(VH,V12,k)
    y=1 ./(1+exp((V12-VH) ./k));
endfunction

//Ca,t+Kir+K,p+L-model
function [Hdot]=HH(t,x,pa)
    Hdot=zeros(4,1);
    Hdot(1)=(1/pa(22))*(-pa(1)*x(2)*x(3)*(x(1)-pa(5)) - pa(2)*xinf(x(1),pa(10),pa(14))*(x(1)-pa(6)) - pa(3)*x(4)*(x(1)-pa(6)) - pa(4)*(x(1)-pa(7)) + I)
    Hdot(2)=(xinf(x(1),pa(8),pa(12))-x(2))/pa(16)
    Hdot(3)=(xinf(x(1),pa(9),pa(13))-x(3))/pa(17)
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
dev=[]
dev1=sigma(A(500:$,1));
dev=[dev1]
for i=2:length(stim)
    dev=[dev sigma(A(5000:$,i))]
end

%ODEOPTIONS=[1,0,0,%inf,0,2,20000,12,5,0,-1,-1];

//Cost function voltage
function y=fct11(pa)
    tmp=0;
    condini = [-53; pa(19); pa(20); pa(21)]
    for i=1:length(stim)
        c=0;
        I=stim(i);
        x=ode(condini,t0,t,HH); 
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-A(k,i))*(V(k)-A(k,i));
        end
        c=sqrt(c/length(t))/dev(i);
        tmp=tmp+c;
    end
    y=tmp/length(stim);
endfunction

///////////////////////////////////////////////////////////////////////
///////////////    Steady-state current cost function    //////////////
///////////////////////////////////////////////////////////////////////

function y=WSS(pa)
    e=0;
    for i=1:length(vecV)
        tmp=0;
        tmp=(Inf(i)-(pa(1)*xinf(vecV(i),pa(8),pa(12))*xinf(vecV(i),pa(9),pa(13))*(vecV(i)-pa(5)) + pa(2)*xinf(vecV(i),pa(10),pa(14))*(vecV(i)-pa(6)) + pa(3)*xinf(vecV(i),pa(11),pa(15))*(vecV(i)-pa(6)) + pa(4)*(vecV(i)-pa(7))))^2
        tmp=tmp/InfSD(i)
        e=e+tmp;
    end
    y=e/length(vecV) 
endfunction

///////////////////////////////////////////////////////////////////////
///////////////    Sum weighted cost function    //////////////////////
///////////////////////////////////////////////////////////////////////

//Weight
wv=100;
wSS=1;

function y=fct(pa)
    y = wv*fct11(pa) + wSS*WSS(pa)
endfunction


//////////////////////////////////////////////
/////////    Parameter estimation    /////////
//////////////////////////////////////////////

function [bM, valBest]=simulation(NP,itermax,F,CR)
    
    D=22;
//    costVec=zeros(1,itermax);
    pop=zeros(D,NP);

    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[0.0001 0.0001 0.0001 0.0001 20  -100 -90 -90 -90 -90 -90 1  -30 -30 1  0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001];
    Xmax=[50     50     50     50     150 -2   30  -2  -2  -2  -2  30 -1  -1  30 20     20     20     0.999 0.999 0.999 10];
    
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
    
    disp(val);
    
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
//    costVec(1)=val(bestIndex);
    
    ////////////////////////
    //// Étape suivante ////
    ////////////////////////
     
    iter=1; // nombre d'itération
    U=zeros(D,NP); // Vecteur intermédiaire perturbé (mutation + crossover)
    tempval=0;
    while iter<itermax
        for j=1:NP
            // ======= Construction de la matrice U = variation différentielle + crossover =======

            // ========= Tirage de 3 individus distincts r1, r2 et r3 suivant le mécanisme de Buhry ========
            if max(val)==min(val) then
                proba=ones(NP,1);
            else
                proba=exp(-5*(val-min(val))/(max(val)-min(val)));//donne un vecteur
                proba=proba/sum(proba);
            end

            r1=j; r2=j; r3=j;
            while (r1==r2 | r1==r3 | r2==r3 | r1==j | r2==j | r3==j)
                r1=find(1==grand(1,'mul',1,proba(1:$-1)))
                r2=find(1==grand(1,'mul',1,proba(1:$-1)))
                r3=find(1==grand(1,'mul',1,proba(1:$-1)))
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
    valBest=val(bestIndex);
    
    // Sauvegarde du meilleur individu
    bM = [];
    bM = pop(:,bestIndex);
    
    disp(val);
endfunction
