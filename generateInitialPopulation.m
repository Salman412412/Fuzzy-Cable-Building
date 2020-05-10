function InitialPopulation = generateInitialPopulation(Population, lb, ub)
%generateInitialPopulation Generates Initial Population
%   This function generates a matrix of number of population row and
%   columns of genes
InitialPopulation = zeros(Population, length(lb));
for pop = 1:Population
    for i = 1:length(lb)
        InitialPopulation(pop,i) = (ub(i)-lb(i))*rand() + lb(i);
    end
end
end
