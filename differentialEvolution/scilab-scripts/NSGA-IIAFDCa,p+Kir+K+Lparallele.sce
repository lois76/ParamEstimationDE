//////////////////////////////////////////////////////////
///////////////     Experimental data      ///////////////
//////////////////////////////////////////////////////////

// Voltage 
a = read("/scilab-scripts/Fig 1A_AFD Current-Clamp Trace.txt",-1,12);
A=a(2489:14988,2:$)*1000;
//a = read("/home/loisse/Documents/FichierScilab/EstimationAFD/Fig1A_AFDCurrentClampTrace.txt",-1,12);
//A=a(2489:14988,2:$)*1000;
t=linspace(0,50,12500);
t0=0;
stim=[-15:5:35];

// Steady-state currents
vecV=[-110:10:50]
Inf=[-68.6 -49.5 -18.2 -5.06 2.19 3.37 2.52 2.68 5.97 14.6 33.4 60.2 85 114 152 208 254]

////////////////////////////////////////////////////////
///////////////    Objective functions    //////////////
////////////////////////////////////////////////////////

function u=xinf(VH,V12,k)
    u = 1 ./(1+exp((V12-VH) ./k));
endfunction

function [Hdot]=HH(t,x,pa)
    Hdot=zeros(4,1);
    Hdot(1)=(1/pa(22))*(-pa(1)*x(2)*(x(1)-pa(5)) - pa(2)*xinf(x(1),pa(9),pa(13))*(x(1)-pa(6)) - pa(3)*x(3)*x(4)*(x(1)-pa(6)) - pa(4)*(x(1)-pa(7)) + I)
    Hdot(2)=(xinf(x(1),pa(8),pa(12))-x(2))/pa(16)
    Hdot(3)=(xinf(x(1),pa(10),pa(14))-x(3))/pa(17)
    Hdot(4)=(xinf(x(1),pa(11),pa(15))-x(4))/pa(18)
endfunction

// Objective function
function f=costFct(pa)
    
    // cost function f1
    c=0;
    condini = [-78; pa(19); pa(20); pa(21)]
    for i=1:11
        I=stim(i);
        x=ode(condini,t0,t,HH); 
        V=x(1,:);
        for k=1:length(t)
            c=c+(V(k)-A(k,i))*(V(k)-A(k,i))
        end
    end
    f1=c/length(t);
    
    // cost function f2
    e=0;
    for i=1:length(vecV)
        e=e+(Inf(i)-(pa(1).*xinf(vecV(i),pa(8),pa(12)).*(vecV(i)-pa(5)) + pa(2).*xinf(vecV(i),pa(9),pa(13)).*(vecV(i)-pa(6)) + pa(3).*xinf(vecV(i),pa(10),pa(14)).*xinf(vecV(i),pa(11),pa(15)).*(vecV(i)-pa(6)) + pa(4).*(vecV(i)-pa(7))))^2
    end
    f2=e/length(vecV) 

    f=[f1, f2];
endfunction

// Borne sup et inf

Xmin=[0.0001 0.0001 0.0001 0.0001 20  -100 -90 -90 -90 -90 -90 1  -30 1  -30 0.0001 0.0001 0.0001 0.001 0.001 0.001 0.001]';
Xmax=[50     50     50     50     150 -2   30  -2  -2  -2  -2  30 -1  30 -1  15     15     15     0.999 0.999 0.999 10]';

// Problem dimension
dim = 22;

// Parameter of the genetic algorithm
funcname='costFct';
PopSize=500;
Proba_cross=0.9;
Proba_mut=0.045 // (=1/22=1/NdVs)
NbGen=5000;
NbCouples=110;
Log=%T;
pressure=0.1;

// Setting parameters of optim_nsga2 function 
ga_params = init_param();

// Parameters to adapt to the shape of the optimization problem
ga_params=add_param(ga_params,'minbound',Xmin);
ga_params=add_param(ga_params,'maxbound',Xmax);
ga_params=add_param(ga_params,'dimension',dim);
ga_params=add_param(ga_params,'beta',0);
ga_params=add_param(ga_params,'delta',0.1);

// Scilab ( http://www.scilab.org/ ) - This file is part of Scilab
// Copyright (C) 2008 - Yann COLLETTE <yann.collette@renault.com>
//
// Copyright (C) 2012 - 2016 - Scilab Enterprises
//
// This file is hereby licensed under the terms of the GNU GPL v2.0,
// pursuant to article 5.3.4 of the CeCILL v.2.1.
// This file was originally licensed under the terms of the CeCILL v2.1,
// and continues to be available under such terms.
// For more information, see the COPYING file which you should have received
// along with this program.

//function [pop_opt, fobj_pop_opt, pop_init, fobj_pop_init, pop_1000, fobj_pop_1000, pop_2000, fobj_pop_2000, pop_3000, fobj_pop_3000, pop_4000, fobj_pop_4000] = optim_nsga2(ga_f, pop_size, nb_generation, p_mut, p_cross, Log, param)

function [pop_opt, fobj_pop_opt, pop_init, fobj_pop_init, pop_1000, fobj_pop_1000, pop_2500, fobj_pop_2500, pop_4000, fobj_pop_4000, pop_4300, fobj_pop_4300, pop_4600, fobj_pop_4600, pop_4800, fobj_pop_4800] = simulation(ga_f, pop_size, nb_generation, p_mut, p_cross, Log, param)

    // isdef vérifie l'existence de la variable (ici en local). Donc ~isdef vérifie la non existence car ~ désigne la négation
    if ~isdef("param", "local") then
        param = [];
    end

    [codage_func, err] = get_param(param, "codage_func", coding_ga_identity);
    [init_func, err] = get_param(param, "init_func", init_ga_default);
    [crossover_func, err] = get_param(param, "crossover_func", crossover_ga_default);
    [mutation_func, err] = get_param(param, "mutation_func", mutation_ga_default);
    [nb_couples, err] = get_param(param, "nb_couples", 100);
    [output_func, err] = get_param(param, 'output_func', output_nsga2_default);

    // si ga_f n'existe pas, alors renvoyé un message d'erreur car celle-ci est obligatoire
    if ~isdef("ga_f", "local") then
        error(gettext("optim_nsga2: ga_f is mandatory"));
    else
        if typeof(ga_f)=="list" then
            deff("y=_ga_f(x)", "y=ga_f(1)(x, ga_f(2:$))");
        else
            deff("y=_ga_f(x)", "y=ga_f(x)");
        end
    end

    // Initialisation de la population
    Pop = init_func(pop_size, param);
    disp(Pop);
    tmp1=Pop;
    pop_init = tmp1;

    // Code the individuals
    Pop = codage_func(Pop, "code", param);
    
    for i=1:length(Pop)
        FObj_Pop(i, :) = _ga_f(Pop(i));
    end
    tmp2=FObj_Pop;
    fobj_pop_init = tmp2;

    // Compute the domination rank
    Rank=DominationRank(FObj_Pop);

    // Compute the crowding distance
    MO_FObj_Pop = FObj_Pop;
    Index    = 1:size(MO_FObj_Pop, 1);
    Crowdist = zeros(size(MO_FObj_Pop, 1), 1);
    for i=1:size(FObj_Pop, 2)
        [tmp, Index_List] = gsort(MO_FObj_Pop(:, i));
        MO_FObj_Pop       = MO_FObj_Pop(Index_List, :);
        Index             = Index(Index_List);
        Crowdist(Index_List(1)) = %inf;
        Crowdist(Index_List($)) = %inf;
        _Max = max(MO_FObj_Pop(:, i));
        _Min = min(MO_FObj_Pop(:, i));
        for j=2:size(MO_FObj_Pop, 1)-1
            Crowdist(Index(j)) = Crowdist(Index(j)) - (MO_FObj_Pop(j+1, i) - MO_FObj_Pop(j-1, i)) / (_Max - _Min);
        end
    end

    // The genetic algorithm
    for It=1:nb_generation
        //
        // Selection
        //
        Indiv1 = list();
        Indiv2 = list();
        for j=1:nb_couples
            // Selection of 2 individuals via binary tournament selection to fill Indiv1
            Index1 = ceil((size(FObj_Pop, 1) - 1)*grand(1, 1, "def")+1);
            Index2 = ceil((size(FObj_Pop, 1) - 1)*grand(1, 1, "def")+1);
            if (Rank(Index1)<Rank(Index2)) | ((Rank(Index1)==Rank(Index2)) & (Crowdist(Index1)>Crowdist(Index2))) then
                Indiv1(j)        = Pop(Index1);
                FObj_Indiv1(j, :) = MO_FObj_Pop(Index1, :);
            else
                Indiv1(j)        = Pop(Index2);
                FObj_Indiv1(j, :) = MO_FObj_Pop(Index2, :);
            end
            // Selection of 2 individuals via binary tournament selection to fill Indiv2
            Index1 = ceil((size(FObj_Pop, 1) - 1)*grand(1, 1, "def")+1);
            Index2 = ceil((size(FObj_Pop, 1) - 1)*grand(1, 1, "def")+1);
            if (Rank(Index1)<Rank(Index2)) | ((Rank(Index1)==Rank(Index2)) & (Crowdist(Index1)>Crowdist(Index2))) then
                Indiv2(j)        = Pop(Index1);
                FObj_Indiv2(j, :) = MO_FObj_Pop(Index1, :);
            else
                Indiv2(j)        = Pop(Index2);
                FObj_Indiv2(j, :) = MO_FObj_Pop(Index2, :);
            end
        end
        //
        // Crossover
        //
        for j=1:nb_couples
            if (p_cross>grand(1, 1, "def")) then
                [x1, x2] = crossover_func(Indiv1(j), Indiv2(j), param);
                Indiv1(j) = x1;
                Indiv2(j) = x2;
                ToCompute_I1(j) = %T;
                ToCompute_I2(j) = %T;
            else
                ToCompute_I1(j) = %F;
                ToCompute_I2(j) = %F;
            end
        end
        //
        // Mutation
        //
        for j=1:nb_couples
            if (p_mut>grand(1, 1, "def")) then
                x1 = mutation_func(Indiv1(j), param);
                Indiv1(j) = x1;
                ToCompute_I1(j) = %T;
            end
            if (p_mut>grand(1, 1, "def")) then
                x2 = mutation_func(Indiv2(j), param);
                Indiv2(j) = x2;
                ToCompute_I2(j) = %T;
            end
        end
        //
        // Computation of the objective functions
        //
        for j=1:length(Indiv1)
            if ToCompute_I1(j) then FObj_Indiv1(j, :) = _ga_f(Indiv1(j)); end
            if ToCompute_I2(j) then FObj_Indiv2(j, :) = _ga_f(Indiv2(j)); end
        end

        // Reinit ToCompute lists
        ToCompute_I1 = ToCompute_I1 & %F;
        ToCompute_I2 = ToCompute_I2 & %F;

        // We merge all the individuals in one list ...
        All_Pop  = lstcat(Pop, Indiv1, Indiv2);
        All_FObj = [FObj_Pop' FObj_Indiv1' FObj_Indiv2']';

        // Compute the domination rank on all the population
        Rank=DominationRank(All_FObj);

        // Compute the crowding distance
        MO_All_FObj = All_FObj;
        Index    = 1:size(MO_All_FObj, 1);
        Crowdist = zeros(size(MO_All_FObj, 1), 1);
        for k=1:size(MO_All_FObj, 2)
            [tmp, Index_List] = gsort(MO_All_FObj(:, k));
            MO_All_FObj = MO_All_FObj(Index_List, :);
            Index = Index(Index_List);
            Crowdist(Index_List(1)) = %inf;
            Crowdist(Index_List($)) = %inf;
            _Max = max(MO_All_FObj(:, k));
            _Min = min(MO_All_FObj(:, k));
            for j=2:size(MO_All_FObj, 1)-1
                Crowdist(Index(j)) = Crowdist(Index(j)) - (MO_All_FObj(j+1, k) - MO_All_FObj(j-1, k)) / (_Max - _Min);
            end
        end
        //
        // Recombination
        //
        // We rank all the individual wrt to the partial order
        for k=1:size(All_FObj, 1)-1
            for j=k+1:size(All_FObj, 1)
                if (Rank(j)<Rank(k)) | ((Rank(j)==Rank(k)) & (Crowdist(j)>Crowdist(k))) then
                    tmp           = Rank(k);
                    Rank(k)       = Rank(j);
                    Rank(j)       = tmp;
                    tmp           = Crowdist(k);
                    Crowdist(k)   = Crowdist(j);
                    Crowdist(j)   = tmp;
                    tmp           = All_Pop(k);
                    All_Pop(k)    = All_Pop(j);
                    All_Pop(j)    = tmp;
                    tmp           = All_FObj(k, :);
                    All_FObj(k, :) = All_FObj(j, :);
                    All_FObj(j, :) = tmp;
                end
            end
        end
        // Extraction and selection of the phenotype
        FObj_Pop = All_FObj(1:pop_size, :);
        // Extraction and selection of the genotype
        Pop = list(All_Pop(1:pop_size));
        // Extraction of the ranks and Crow distance
        Rank     = Rank(1:pop_size);
        Crowdist = Crowdist(1:pop_size);
        
        // Save population and associated objective fonctions
        if It==1000 then
            disp(Pop);
            npop=length(Pop);
            pop_1000=matrix(list2vec(Pop),dim,npop);
            fobj_pop_1000=FObj_Pop;
        end
        if It==2500 then
            disp(Pop);
            npop=length(Pop);
            pop_2500=matrix(list2vec(Pop),dim,npop);
            fobj_pop_2500=FObj_Pop;
        end
        if It==4000 then
            disp(Pop);
            npop=length(Pop);
            pop_4000=matrix(list2vec(Pop),dim,npop);
            fobj_pop_4000=FObj_Pop;
        end
        if It==4300 then
            disp(Pop);
            npop=length(Pop);
            pop_4300=matrix(list2vec(Pop),dim,npop);
            fobj_pop_4300=FObj_Pop;
        end
        if It==4600 then
            disp(Pop);
            npop=length(Pop);
            pop_4600=matrix(list2vec(Pop),dim,npop);
            fobj_pop_4600=FObj_Pop;
        end
        if It==4800 then
            disp(Pop);
            npop=length(Pop);
            pop_4800=matrix(list2vec(Pop),dim,npop);
            fobj_pop_4800=FObj_Pop;
        end
        if It==nb_generation then
            disp(Pop);
            npop=length(Pop);
            pop_opt=matrix(list2vec(Pop),dim,npop);
            fobj_pop_opt=FObj_Pop;
        end

        if (Log) then
            stop = output_func(i, nb_generation, Pop, FObj_Pop, param);
            if stop then
                break
            end
        end
    end
endfunction

//// Performing optimization
//printf("Performing optimization:");
[pop_opt, fobj_pop_opt, pop_init, fobj_pop_init, pop_3, fobj_pop_3, pop_5, fobj_pop_5] = simulation(costFct,PopSize,NbGen,Proba_mut,Proba_cross,Log,ga_params);


