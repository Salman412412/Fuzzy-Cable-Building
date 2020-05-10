clear,clc, close all

%---------------------------------------
% --- Define Building for Simulation ---
%---------------------------------------
No_bld     = 20;			% 20 Story Building

%---------------------------------------------
% --- Build Benchmark Model for Simulation ---
%---------------------------------------------
sprintf(' ---> BUILDING MODEL FOR SIMULATION')
Bld_NLBM

%----------------------------------------------
% --- Simulate El Centro Earthquake Records ---
%----------------------------------------------
OPTIONS = simset('solver','ode5','FixedStep',dt_cal);
isim = 0;
tf       = 100;            % Duration (sec)

% For training we use Elcentro by intensity of 1.5
intensity = 1.5;


% These are the genetic algorithm parametes set by Bozorgvar essay Table 4
fun = @(xsample)costs_2016(xsample,tf,OPTIONS);
opts = optimoptions('ga');
opts.PopulationSize = 16;  % This should be even
opts.EliteCount = 2;
opts.MaxGenerations = 10;
opts.CreationFcn = 'gacreationlinearfeasible';
opts.MigrationDirection = 'both';
opts.FitnessLimit = 0;
opts.MutationFcn = 'mutationadaptfeasible';
CrossOverFactor = 0.8;
MutationFactor = 0.2;


global Params x1;    % Params is what should be trained to minimize J1
% Param 1 to 10 are sigma1, c1, ..., sigma5, c5 of acc
% Param 11 to 20 are sigma1, c1, ..., sigma5, c5 of drift
% Param 21 to 95 are p1, q1, r1, ..., p25, q25, r25
lb = zeros(1,95);   % lb is the minimum of the possible values in Param
ub = zeros(1,95);   % ub is the maximum of the possible values in Param

for i = 1:2:20      % This section says sigmas should be between 0.5 and 2
    lb(i) = 0.5;
    ub(i) = 2;
end
for i = 2:2:20      % This section says cs should be between -10 and 10
    lb(i) = -10;
    ub(i) = 10;
end
for i = 21:95       % This section says p, q, and r should be between 0 and 10
    lb(i) = 0;
    ub(i) = 10;
end

warning('off')

% Genetic Algorithm Initialization
Params = generateInitialPopulation(opts.PopulationSize, lb, ub);
J1 = zeros(1,opts.PopulationSize);
bestJ1 = Inf*ones(1,opts.MaxGenerations);
allbestJ1 = Inf;

tic
for Generation = 1:opts.MaxGenerations
    
    % Calculation of the J1 of the population
    for pop = 1:opts.PopulationSize
        x1 = zeros(tf*100+1, 96);
        x1(:,1) = (0:0.01:tf)';
        for i = 1:(tf*100+1)
            x1(i,2:96) = (Params(pop,:));
        end
        J1(pop) = costs_2016(tf, OPTIONS);
        if (Generation == 1) && (pop == 1)
            time = toc;
            remaining = opts.MaxGenerations * opts.PopulationSize * time;
            days = floor(remaining/(3600*24));
            remaining = mod(remaining, (3600*24));
            hours = floor(remaining/(3600));
            remaining = mod(remaining, (3600));
            minutes = floor(remaining/(60));
            seconds = mod(remaining, (60));
            fprintf('The remaining time (first estimation): %d days, %d hours, %d minutes, and %2.1f seconds\n', days, hours, minutes, seconds);
        end
    end
    
    % Calculating Elites
    [J1Sorted, Indices] = sort(J1);
    Elites = Params(Indices(1:opts.EliteCount),:);
    
    % Check for the best J1 of the all generations
    if J1Sorted(1)<allbestJ1
        allbestJ1 = J1Sorted(1);
        bestGenes = Elites(1,:);
    end
    
    % Printing the results
    hold on
    plot(Generation, sum(J1)/length(J1), '.b');
    plot(Generation, J1Sorted(1), '*r');
    axis([1, opts.MaxGenerations, 0.2, 1.2])
    drawnow;
    
    % Calculates Scale Fitness
    Fitness = zeros(1,opts.PopulationSize);
    n = 1;
    for i = Indices
        Fitness(i) = 1/sqrt(n);
        n = n + 1;
    end
    
    % Finding Parents
    NoParents = opts.PopulationSize - opts.EliteCount;
    Parents = zeros(NoParents,95);
    TotalFitness = sum(Fitness);
    Scale = TotalFitness / NoParents;
    FitnessSorted = Fitness(Indices);
    i = 1; flag = false;
    while i <= NoParents
        MinFitness = i * Scale;
        FitnessSum = 0;
        for pop = opts.PopulationSize:-1:2
            FitnessSum = FitnessSum + FitnessSorted(pop);
            if FitnessSum >= MinFitness
                Parents(i,:) = Params(Indices(pop),:);
                i = i + 1;
                break
            end
        end
        pop = 1;
        Parents(i,:) = Params(Indices(pop),:);
        i = i + 1;
    end
    
    % Mating
    ParentIndices = randperm(NoParents);
    Childs = zeros(NoParents, 95);
    for child = 1:2:NoParents
        randomVector = (rand(1,95)<=CrossOverFactor);
        for gene = 1:95
            if randomVector(gene)
                Childs(child,gene) = Parents(ParentIndices(child),gene);
                Childs(child+1,gene) = Parents(ParentIndices(child+1),gene);
            else
                Childs(child,gene) = Parents(ParentIndices(child+1),gene);
                Childs(child+1,gene) = Parents(ParentIndices(child),gene);
            end
        end
    end
    
    % Collecting next generation
    NextGen = [Elites;Childs];
    for pop = 1:opts.PopulationSize
        MutationVector = (rand(1,95)<=MutationFactor);
        for gene = 1:95
            if MutationVector(gene)
                NextGen(pop,gene) = (ub(gene)-lb(gene))*rand() + lb(gene);
            end
        end
    end
    
    % Replace the genes with the Generation
    Params = NextGen;
    
    % Calculate and show the remaining time
    time = toc;
    remaining = (opts.MaxGenerations-Generation)/Generation*time;
    days = floor(remaining/(3600*24));
    remaining = mod(remaining, (3600*24));
    hours = floor(remaining/(3600));
    remaining = mod(remaining, (3600));
    minutes = floor(remaining/(60));
    seconds = mod(remaining, (60));
    fprintf('The remaining time: %d days, %d hours, %d minutes, and %2.1f seconds\n', days, hours, minutes, seconds);
end
fprintf('The best J1 is: %0.4f\n', allbestJ1);
save('TheBest.mat', 'bestGenes', 'allbestJ1')