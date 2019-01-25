//////////////////////////////////////////////////////////
///////////////     Simulation modèle      ///////////////
//////////////////////////////////////////////////////////

/// Potentiel d'équilibre ///
ENa=120; EK=-12; EL=10;

/// Conductance ///
gNa=3; gK=4; gL=3;


function y=xinf(V,V12,k)
    y=1/(1+exp((V12-V)/k));
endfunction

function [Hdot]=HH(t,x)
    Hdot=zeros(4,1); // définition de la taille
    Hdot(1)=-gNa*x(2)*x(3)*(x(1)-ENa)-gK*x(4)*(x(1)-EK)-gL*(x(1)-EL)+I
    Hdot(2)=(xinf(x(1),-40,15)-x(2))/5
    Hdot(3)=(xinf(x(1),-62,-7)-x(3))/9
    Hdot(4)=(xinf(x(1),-53,11)-x(4))/3
endfunction

/// Stimulations appliquées au neurone ///
stim=[-15:5:35];

/// Construction des solutions ///
t0=0;
t=0:0.05:50;
t=t';
a=zeros(length(t),12);
a(:,1)=t;
for i=1:11
    I=stim(i);
    x=ode([-35;0.1;0.3;0.2],t0,t,HH);
    x1=x(1,:);
    x1=x1';
    a(:,i+1)=x1;
end


///////////////////////////////////////////////////////////////////////////
///////////////     Récupération données expérimentales     ///////////////
///////////////////////////////////////////////////////////////////////////

// Récupération des données expérimentales sous forme de matrice
//A = read("/home/loisse/Documents/FichierScilab/Fig1A_AIYCurrentClampTrace2.txt",-1,12);
//a=A*1000;
//t=a(:,1)/1000;
//t0=t(1);
//stim=[-15:5:35]; // vecteur des différents courants d'injections


/////////////////////////////////////////////////////////////////////////
//////////////////// Définitions fonctions de canaux ////////////////////
/////////////////////////////////////////////////////////////////////////

function y=xinf(V,V12,k)
    y=1/(1+exp((V12-V)/k));
endfunction

function y=ICa1(x1,x2,V,j)
    y=popTemp(3,j)*x1*x2*(V-popTemp(6,j))
endfunction

function y=ICa2(x1,V,j)
    y=popTemp(3,j)*x1*(V-popTemp(6,j))
endfunction

function y=ICa3(V,j,V12,k)
    y=popTemp(3,j)*xinf(V,V12,k)*(V-popTemp(6,j))
endfunction

function y=IK1(x3,x4,V,j)
    y=popTemp(4,j)*x3*x4*(V-popTemp(7,j))
endfunction

function y=IK2(x3,V,j)
    y=popTemp(4,j)*x3*(V-popTemp(7,j))
endfunction

function y=IK3(V,j,V12,k)
    y=popTemp(4,j)*xinf(V,V12,k)*(V-popTemp(7,j))
endfunction

function y=IL(V,j) 
    y=popTemp(5,j)*(V-popTemp(8,j))
endfunction



//////////////////////////////////////////////////////////////////////////////
/////// Définitions systèmes d'équations en fonction du type de canaux ///////
//////////////////////////////////////////////////////////////////////////////

function [Hdot]=HH11(t,x,j)
    Hdot=zeros(5,1);
    Hdot(1)=(1/popTemp(21,j))*(-ICa1(x(2),x(3),x(1),j)-IK1(x(4),x(5),x(1),j)-IL(x(1),j)+I)
    Hdot(2)=(xinf(x(1),popTemp(9,j),popTemp(13,j))-x(2))/popTemp(17,j)
    Hdot(3)=(xinf(x(1),popTemp(10,j),popTemp(14,j))-x(3))/popTemp(18,j)
    Hdot(4)=(xinf(x(1),popTemp(11,j),popTemp(15,j))-x(4))/popTemp(19,j)
    Hdot(5)=(xinf(x(1),popTemp(12,j),popTemp(16,j))-x(5))/popTemp(20,j)
endfunction

function [Hdot]=HH12(t,x,j)
    Hdot=zeros(4,1);
    Hdot(1)=(1/popTemp(21,j))*(-ICa1(x(2),x(3),x(1),j)-IK2(x(4),x(1),j)-IL(x(1),j)+I)
    Hdot(2)=(xinf(x(1),popTemp(9,j),popTemp(13,j))-x(2))/popTemp(17,j)
    Hdot(3)=(xinf(x(1),popTemp(10,j),popTemp(14,j))-x(3))/popTemp(18,j)
    Hdot(4)=(xinf(x(1),popTemp(11,j),popTemp(15,j))-x(4))/popTemp(19,j)
endfunction

function [Hdot]=HH13(t,x,j)
    Hdot=zeros(3,1);
    Hdot(1)=(1/popTemp(21,j))*(-ICa1(x(2),x(3),x(1),j)-IK3(x(1),j,popTemp(11,j),popTemp(15,j))-IL(x(1),j)+I)
    Hdot(2)=(xinf(x(1),popTemp(9,j),popTemp(13,j))-x(2))/popTemp(17,j)
    Hdot(3)=(xinf(x(1),popTemp(10,j),popTemp(14,j))-x(3))/popTemp(18,j)
endfunction

function [Hdot]=HH21(t,x,j)
    Hdot=zeros(4,1);
    Hdot(1)=(1/popTemp(21,j))*(-ICa2(x(2),x(1),j)-IK1(x(3),x(4),x(1),j)-IL(x(1),j)+I)
    Hdot(2)=(xinf(x(1),popTemp(9,j),popTemp(13,j))-x(2))/popTemp(17,j)
    Hdot(3)=(xinf(x(1),popTemp(11,j),popTemp(15,j))-x(3))/popTemp(19,j)
    Hdot(4)=(xinf(x(1),popTemp(12,j),popTemp(16,j))-x(4))/popTemp(20,j)
endfunction

function [Hdot]=HH22(t,x,j)
    Hdot=zeros(3,1);
    Hdot(1)=(1/popTemp(21,j))*(-ICa2(x(2),x(1),j)-IK2(x(2),x(1),j)-IL(x(1),j)+I)
    Hdot(2)=(xinf(x(1),popTemp(9,j),popTemp(13,j))-x(2))/popTemp(17,j)
    Hdot(3)=(xinf(x(1),popTemp(11,j),popTemp(15,j))-x(3))/popTemp(19,j)
endfunction

function [Hdot]=HH23(t,x,j)
    Hdot=zeros(2,1);
    Hdot(1)=(1/popTemp(21,j))*(-ICa2(x(2),x(1),j)-IK3(x(1),j,popTemp(11,j),popTemp(15,j))-IL(x(1),j)+I)
    Hdot(2)=(xinf(x(1),popTemp(9,j),popTemp(13,j))-x(2))/popTemp(17,j)
endfunction

function [Hdot]=HH31(t,x,j)
    Hdot=zeros(3,1);
    Hdot(1)=(1/popTemp(21,j))*(-ICa3(x(1),j,popTemp(9,j),popTemp(13,j))-IK1(x(2),x(3),x(1),j)-IL(x(1),j)+I)
    Hdot(2)=(xinf(x(1),popTemp(11,j),popTemp(15,j))-x(2))/popTemp(19,j)
    Hdot(3)=(xinf(x(1),popTemp(12,j),popTemp(16,j))-x(3))/popTemp(20,j)
endfunction

function [Hdot]=HH32(t,x,j)
    Hdot=zeros(2,1);
    Hdot(1)=(1/popTemp(21,j))*(-ICa3(x(1),j,popTemp(9,j),popTemp(13,j))-IK2(x(2),x(1),j)-IL(x(1),j)+I)
    Hdot(2)=(xinf(x(1),popTemp(11,j),popTemp(15,j))-x(2))/popTemp(19,j)
endfunction

function [Hdot]=HH33(t,x,j)
    Hdot=zeros(1,1);
    Hdot(1)=(1/popTemp(21,j))*(-ICa3(x(1),j,popTemp(9,j),popTemp(13,j))-IK3(x(1),j,popTemp(11,j),popTemp(15,j))-IL(x(1),j)+I)
endfunction


//////////////////////////////////////////////////////////////////////////////
///////// Définitions fonctions coûts associées aux types de canaux  /////////
//////////////////////////////////////////////////////////////////////////////

function y=fct11(W)
    c=0;
    popTemp = W;
    condini = [V0; popTemp(22,j); popTemp(23,j); popTemp(24,j); popTemp(25,j)]
    for i=1:11
        I=stim(i);
        x=ode(condini,t0,t,HH11); // Quelles conditions initiales mettre ?
        V=x(1,:); //length(V) pour voir s'il est de la même taille que a(k,i)
        for k=1:length(t)
            c=c+(V(k)-a(k,i+1))*(V(k)-a(k,i+1))
        end
    end
    y=c;
endfunction

function y=fct12(W)
    c=0;
    popTemp = W;
    condini = [V0; popTemp(22,j); popTemp(23,j); popTemp(24,j)]
    for i=1:11
        I=stim(i)
        x=ode(condini,t0,t,HH12); // Quelles conditions initiales mettre ?
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-a(k,i+1))*(V(k)-a(k,i+1))
        end
    end
    y=c;
endfunction

function y=fct13(W)
    c=0;
    popTemp = W;
    condini = [V0; popTemp(22,j); popTemp(23,j)]
    for i=1:11
        I=stim(i)
        x=ode(condini,t0,t,HH13); // Quelles conditions initiales mettre ?
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-a(k,i+1))*(V(k)-a(k,i+1))
        end
    end
    y=c;
endfunction

function y=fct21(W)
    c=0;
    popTemp = W;
    condini = [V0; popTemp(22,j); popTemp(24,j); popTemp(25,j)]
    for i=1:11
        I=stim(i)
        x=ode(condini,t0,t,HH21); // Quelles conditions initiales mettre ?
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-a(k,i+1))*(V(k)-a(k,i+1))
        end
    end
    y=c;
endfunction

function y=fct22(W)
    c=0;
    popTemp = W;
    condini = [V0; popTemp(22,j); popTemp(24,j)]
    for i=1:11
        I=stim(i)
        x=ode(condini,t0,t,HH22); // Quelles conditions initiales mettre ?
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-a(k,i+1))*(V(k)-a(k,i+1))
        end
    end
    y=c;
endfunction

function y=fct23(W)
    c=0;
    popTemp = W;
    condini = [V0; popTemp(22,j)]
    for i=1:11
        I=stim(i)
        x=ode(condini,t0,t,HH23); // Quelles conditions initiales mettre ?
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-a(k,i+1))*(V(k)-a(k,i+1))
        end
    end
    y=c;
endfunction

function y=fct31(W)
    c=0;
    popTemp = W;
    condini = [V0; popTemp(24,j); popTemp(25,j)]
    for i=1:11
        I=stim(i)
        x=ode(condini,t0,t,HH31); // Quelles conditions initiales mettre ?
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-a(k,i+1))*(V(k)-a(k,i+1))
        end
    end
    y=c;
endfunction

function y=fct32(W)
    c=0;
    popTemp = W;
    condini = [V0; popTemp(24,j)]
    for i=1:11
        I=stim(i)
        x=ode(condini,t0,t,HH32); // Quelles conditions initiales mettre ?
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-a(k,i+1))*(V(k)-a(k,i+1))
        end
    end
    y=c;
endfunction

function y=fct33(W)
    c=0;
    popTemp = W;
    condini = V0;
    for i=1:11
        I=stim(i)
        x=ode(condini,t0,t,HH33); // Quelles conditions initiales mettre ?
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-a(k,i+1))*(V(k)-a(k,i+1))
        end
    end
    y=c;
endfunction


function main(NP,itermax,F,CR)

    t1=getdate();
    V0=-35; // potentiel de repos
    D=25; // taille des individus
    costVec=zeros(1,itermax); // sauvegarde le coût du meilleur individu à chaque itération
    
    ///////////////////////////////////////////////////////
    //// Vecteurs de contraintes borne minimum/maximum ////
    ///////////////////////////////////////////////////////

    Xmin=[1 1 1 1 1 100 -50 0 -77 -77 -77 -77 1 -15 1 -15 1 1 1 1 1];
    Xmax=[3 3 10 10 10 140 0 20 -1 -1 -1 -1 20 -1 20 -1 20 20 20 20 10];
    //Xmin=[1 1 1 1 1 50 -100 -80 -77 -77 -77 -77 1 -15 1 -15 1 1 1 1 1];
    //Xmax=[3 3 10 10 10 140 -50 -40 -1 -1 -1 -1 20 -1 20 -1 20 20 20 20 15];

    /////////////////////////////////////////
    //// Initialisation de ma population ////
    /////////////////////////////////////////

    pop=zeros(D,NP);
    for j=1:NP
//        pop(1,j)=floor(Xmin(1)+Xmax(1)*rand());
//        pop(2,j)=floor(Xmin(2)+Xmax(2)*rand());
        pop(1,j)=2;
        pop(2,j)=1;
        for i=3:(D-4)
            pop(i,j)=floor(Xmin(i)+(Xmax(i)-Xmin(i))*rand());
        end
        for i=(D-3):D
            pop(i,j)=rand();
        end
    end

    popInit=pop; // matrice population initiale récupérée pour l'écrire dans un fichier .txt
    disp(popInit);

    //////////////////////////////////////////////////////////////
    //// Évaluation du meilleur individu après initialisation ////
    //////////////////////////////////////////////////////////////

    val = zeros(NP,1); // tableau avec le coût de chacun des individus

    // ---------- Evalutation du premier individu après initialisation ----------
    j=1;
    if pop(1,j)==1 & pop(2,j)==1 then
        val(j)=fct11(pop);

    elseif pop(1,j)==1 & pop(2,j)==2 then
        val(j)=fct12(pop);

    elseif pop(1,j)==1 & pop(2,j)==3 then
        val(j)=fct13(pop);

    elseif pop(1,j)==2 & pop(2,j)==1 then
        val(j)=fct21(pop);

    elseif pop(1,j)==2 & pop(2,j)==2 then
        val(j)=fct22(pop);

    elseif pop(1,j)==2 & pop(2,j)==3 then
        val(j)=fct23(pop);

    elseif pop(1,j)==3 & pop(2,j)==1 then
        val(j)=fct31(pop);

    elseif pop(1,j)==3 & pop(2,j)==2 then
        val(j)=fct32(pop);

    elseif pop(1,j)==3 & pop(2,j)==3 then
        val(j)=fct33(pop);
    end//ferme le if

    // ------------- Évaluation des membres restants -------------

    for j=2:NP
        if pop(1,j)==1 & pop(2,j)==1 then
            val(j)=fct11(pop);

        elseif pop(1,j)==1 & pop(2,j)==2 then
            val(j)=fct12(pop);

        elseif pop(1,j)==1 & pop(2,j)==3 then
            val(j)=fct13(pop);

        elseif pop(1,j)==2 & pop(2,j)==1 then
            val(j)=fct21(pop);

        elseif pop(1,j)==2 & pop(2,j)==2 then
            val(j)=fct22(pop);

        elseif pop(1,j)==2 & pop(2,j)==3 then
            val(j)=fct23(pop);

        elseif pop(1,j)==3 & pop(2,j)==1 then
            val(j)=fct31(pop);

        elseif pop(1,j)==3 & pop(2,j)==2 then
            val(j)=fct32(pop);

        elseif pop(1,j)==3 & pop(2,j)==3 then
            val(j)=fct33(pop);
        end//ferme le if
    end //ferme le for j=2:NP
    
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
    costVec(1)=val(bestIndex);

    //// =============== ÉTAPE SUIVANTE ================
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
            for r=1:(D-4)
                V(r)=floor(V(r));
            end
            // ----- Condition pour que tous les termes soient différents de 0 -----
//            if V(1)<1 then V(1)=1;
//            elseif V(1)>3 then V(1)=3;
//            end
//            if V(2)<1 then V(2)=1;
//            elseif V(2)>3 then V(2)=3;
//            end
            V(1)=2;
            V(2)=1;
            if V(3)<=1 then V(3)=1; end
            if V(4)<=1 then V(4)=1; end
            if V(5)<=1 then V(5)=1; end
            if V(3)>10 then V(3)=10; end
            if V(4)>10 then V(4)=10; end
            if V(5)>10 then V(5)=10; end
            if V(6)<100 then V(6)=100; end
            if V(6)>140 then V(6)=140; end
            //        if V(7)<-100 then V(7)=-100; end
            //        if V(7)>-50 then V(7)=-50; end
            if V(7)<-50 then V(7)=-50; end
            if V(7)>-1 then V(7)=-1; end
            //        if V(8)<-80 then V(8)=-80; end
            //        if V(8)>-40 then V(8)=-40; end
            if V(8)<1 then V(8)=1; end
            if V(8)>20 then V(7)=20; end
            if -1 <= V(9) then V(9)=-1; end
            if -1 <= V(10) then V(10)=-1; end
            if -1 <= V(11) then V(11)=-1; end
            if -1 <= V(12) then V(12)=-1; end
            if V(9)<-77 then V(9)=-77; end
            if V(10)<-77 then V(10)=-77; end
            if V(11)<-77 then V(11)=-77; end
            if V(12)<-77 then V(12)=-77; end
            if V(13) <= 1 then V(13)=1; end
            if V(13)>20 then V(13)=20;end
            if -1 <= V(14) then V(14)=-1; end
            if V(14)<-15 then V(14)=-15;end
            if V(15)<=1 then V(15)=1; end
            if V(15)>20 then V(15)=20;end
            if -1<=V(16) then V(16)=-1; end
            if V(16)<-15 then V(16)=-15;end
            if V(17)<=1 then V(17)=1; end
            if V(17)>20 then V(17)=20; end
            if V(18)<=1 then V(18)=1; end
            if V(18)>20 then V(18)=20; end
            if V(19)<=1 then V(19)=1; end
            if V(19)>20 then V(19)=20; end
            if V(20)<=1 then V(20)=1; end
            if V(20)>20 then V(20)=20; end
            if V(21)<=1 then V(21)=1; end
            if V(21)>15 then V(21)=15; end
            if V(22)>0.9 then V(22)=0.9; end
            if V(22)<0.1 then V(22)=0.1; end
            if V(23)>0.9 then V(23)=0.9; end
            if V(23)<0.1 then V(23)=0.1; end
            if V(24)>0.9 then V(24)=0.9; end
            if V(24)<0.1 then V(24)=0.1; end
            if V(25)>0.9 then V(25)=0.9; end
            if V(25)<0.1 then V(25)=0.1; end
            // ======== Crossover ========
            for i=1:D
                if rand()<CR then
                    U(i,j)=V(i);
                else
                    U(i,j)=pop(i,j);
                end
            end
        end // fin for j=1:NP
        // ====== Fin de construction de la matrice U =====

        // ======== Sélection ========
        // !!!!! ATTENTION À LA NOTATION j QUI EST DÉJÀ PRISE JUSTE AU DESSUS ! OU PEUT-ÊTRE QUE C'EST VRAIMENT LA BONNE CAR j EST GLOBAL ET SERA PRISE PAR LA FONCTION fct11, fct12, ...
        for j=1:NP
            if U(1,j)==1 & U(2,j)==1 then
                tempval=fct11(U);

            elseif pop(1,j)==1 & pop(2,j)==2 then
                tempval=fct12(U);

            elseif pop(1,j)==1 & pop(2,j)==3 then
                tempval=fct13(U);

            elseif pop(1,j)==2 & pop(2,j)==1 then
                tempval=fct21(U);

            elseif pop(1,j)==2 & pop(2,j)==2 then
                tempval=fct22(U);

            elseif pop(1,j)==2 & pop(2,j)==3 then
                tempval=fct23(U);

            elseif pop(1,j)==3 & pop(2,j)==1 then
                tempval=fct31(U);

            elseif pop(1,j)==3 & pop(2,j)==2 then
                tempval=fct32(U);

            elseif pop(1,j)==3 & pop(2,j)==3 then
                tempval=fct33(U);
            end
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
    disp(pop);
    disp(val);
    
    // Détermination des compteurs des types de canaux
    compt11=0; compt12=0; compt13=0;
    compt21=0; compt22=0; compt23=0;
    compt31=0; compt32=0; compt33=0;
    for j=1:NP
        if     pop(1,j)==1 & pop(2,j)==1 then compt11=compt11+1;
        elseif pop(1,j)==1 & pop(2,j)==2 then compt12=compt12+1;
        elseif pop(1,j)==1 & pop(2,j)==3 then compt13=compt13+1;
        elseif pop(1,j)==2 & pop(2,j)==1 then compt21=compt21+1;
        elseif pop(1,j)==2 & pop(2,j)==2 then compt22=compt22+1;
        elseif pop(1,j)==2 & pop(2,j)==3 then compt23=compt23+1;
        elseif pop(1,j)==3 & pop(2,j)==1 then compt31=compt31+1;
        elseif pop(1,j)==3 & pop(2,j)==2 then compt32=compt32+1;
        elseif pop(1,j)==3 & pop(2,j)==3 then compt33=compt33+1;
        end
    end
    compt=[];
    compt(1)=compt11; compt(2)=compt12; compt(3)=compt13;
    compt(4)=compt21; compt(5)=compt22; compt(6)=compt23;
    compt(7)=compt31; compt(8)=compt32; compt(9)=compt33;
    disp(compt);
    
    // Détermination de l'inidice du meilleur individu
    bestIndex=1;
    for b=2:NP
        if val(b)<val(bestIndex) then bestIndex=b; end
    end
    disp(bestIndex);
    disp(val(bestIndex));
    
    // Sauvegarde du meilleur individu
    bestMember = [];
    bestMember = pop(:,bestIndex);
    disp(bestMember);
    
    t2=getdate();
    E=etime(t2,t1);
    
    iterVec=1:1:itermax
    plot(iterVec,costVec,2)
//    save(NP,itermax,F,CR,compt,val,pop,popInit,bestIndex,bestMember,E)
endfunction


//function save(NP,itermax,F,CR,compt,val,pop,popInit,bestIndex,bestMember,E)
//    mkdir('/home/naudinl/Documents/DonneesDE/' + string(NP) + "_" + string(itermax) + "_" + string(F) + "_" + string(CR));
//    write('/home/naudinl/Documents/DonneesDE/' + string(NP) + "_" + string(itermax) + "_" + string(F) + "_" + string(CR) + '/popInit.txt',popInit);
//    write('/home/naudinl/Documents/DonneesDE/' + string(NP) + "_" + string(itermax) + "_" + string(F) + "_" + string(CR) + '/popFinal.txt',pop);
//    write('/home/naudinl/Documents/DonneesDE/' + string(NP) + "_" + string(itermax) + "_" + string(F) + "_" + string(CR) + '/cost.txt',val);
//    write('/home/naudinl/Documents/DonneesDE/' + string(NP) + "_" + string(itermax) + "_" + string(F) + "_" + string(CR) + '/countType.txt',compt,"(I4)");
//    write('/home/naudinl/Documents/DonneesDE/' + string(NP) + "_" + string(itermax) + "_" + string(F) + "_" + string(CR) + '/bestIndex.txt',bestIndex,"(I4)");
//    write('/home/naudinl/Documents/DonneesDE/' + string(NP) + "_" + string(itermax) + "_" + string(F) + "_" + string(CR) + '/bestMember.txt',bestMember);
//    write('/home/naudinl/Documents/DonneesDE/' + string(NP) + "_" + string(itermax) + "_" + string(F) + "_" + string(CR) + '/time.txt',E)
//endfunction

main(60,150,0.5,0.9)

//taillesPop = [200 300]
//iterations = [100 300]
//mutations = [0.2 0.5 1]
//CRs = [0.75 0.85 0.95]
////
//taillesPop = [10 20] // 10;
//iterations = [10 20] //15;
//mutations = [0.5] //0.5;
//CRs = [0.85 0.9] 
//
//
//t3=getdate();
//parallel_run(taillesPop,iterations,mutations,CRs,main)
//t4=getdate();
//E=etime(t4,t3);


//S=zeros(4,length(taillesPop)*length(iterations)*length(mutations)*length(CRs))
//j=1;
//for NP = taillesPop
//    for itermax = iterations
//       for F = mutations
//           for CR = CRs
////               disp(j)
////               S(1,j) = NP
////               S(2,j) = itermax
////               S(3,j) = F
////               S(4,j) = CR  
////               j=j+1;
//                main(NP,itermax,F,CR)
//           end
//       end
//    end
//end





///////////////////////////////////////////////////////////////////////
///////////////////////// Analyse des données /////////////////////////
///////////////////////////////////////////////////////////////////////

//mkdir('/home/naudinl/Documents/DonneesDE/DataAnalysis');
//write('/home/naudinl/Documents/DonneesDE/countTypeMean.txt',);
//countTypeMean=zeros(9,1);
//for NP = taillesPop
//    for itermax = iterations
//       for F = mutations
//           for CR = CRs
//                ctm=read('/home/naudinl/Documents/DonneesDE/' + string(NP) + "_" + string(itermax) + "_" + string(F) + "_" + string(CR) + '/countType.txt',9,1);
//                countTypeMean=countTypeMean+ctm;
//           end
//        end
//    end
//end



