function NSGAIIDFCS(Global)
% <algorithm> <N>
% Nondominated sorting genetic algorithm II

%------------------------------- Reference --------------------------------
% K. Deb, A. Pratap, S. Agarwal, and T. Meyarivan, A fast and elitist
% multiobjective genetic algorithm: NSGA-II, IEEE Transactions on
% Evolutionary Computation, 2002, 6(2): 182-197.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2018-2019 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Generate random population
    Population = Global.Initialization();
%     [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);
    [FrontNo,~] = NDSort(Population.objs,Population.cons,Global.N);  %%∑«÷ß≈‰≈≈–Ú
    CrowdDis = CrowdingDistance(Population.objs,FrontNo);     %%”µº∑æ‡¿Î
        Zmin         = min(Population(all(Population.cons<=0,2)).objs,[],1);

    %% Optimization
    while Global.NotTermination(Population)
        MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        Offspring  = GA(Population(MatingPool));
        
        Zmin       = min([Zmin;Offspring(all(Offspring.cons<=0,2)).objs],[],1);
        Population = EnvironmentalSelectionDFCS([Population,Offspring],Global.N,Zmin); 
%         [Population,FrontNo,CrowdDis] = EnvironmentalSelection([Population,Offspring],Global.N);
        [FrontNo,~] = NDSort(Population.objs,Population.cons,Global.N);  %%∑«÷ß≈‰≈≈–Ú
        CrowdDis = CrowdingDistance(Population.objs,FrontNo);     %%”µº∑æ‡¿Î
    end
end