import random # Used to generate random numbers

def init_population (pop_size : int, chrom_len : int, bins : int) -> list:
    '''
    This function initialises a population of solutions to a BPP, using the given parameters.

    :pop_size: The size of the population (number of solutions) to initialise.
    :chrom_len: The length of each solution chromosome. Equal to the number of items in the BPP.
    :bins: The number of bins in the BPP. Defines the domain of the chromosome genes.

    Returns a list containing the population of solutions.
    '''
    population = list()

    for x in range(pop_size):
        # The outer loop iterates once for each solution to be added to the population
        chromosome = []

        for y in range(chrom_len):
            # The inner loop iterates once for each gene/item in a solution chromosome
            chromosome.append(random.randint(1, bins)) # Each gene is a random integer between 1 and the number of bins
        
        population.append(chromosome) # Adds the solution chromosome to the population

    return population


def calc_fitness (chromsome : list, bins : int, which_bpp : bool) -> float:
    '''
    This function calculates the fitness of a given solution.

    :chromosome: The solution chromosome to calculate the fitness of.
    :bins: The number of bins in the BPP.
    :which_bpp: Set to True for coursework BPP1, or False for BPP2.

    Returns an integer denoting the fitness of the solution.
    '''
    bin_weights = [0] * bins # This list will hold the weight of each bin
    
    for gene_index in range(len(chromsome)):
        # This loop adds the weight of each item to the correct bin
        if which_bpp == True: 
            bin_weights[chromsome[gene_index] - 1] += gene_index + 1 # For BPP1
        else:
            bin_weights[chromsome[gene_index] - 1] += ((gene_index + 1) ** 2) / 2 # For BPP2

    fitness = 100 / (1 + (max(bin_weights) - min(bin_weights)))

    return fitness


def select_parents (pop_fitness : list, tour_size : int) -> tuple:
    '''
    This function uses tournament selection in order to select two parent solutions from the population.

    :pop_fitness: The list of the population's fitness values.
    :tour_size: The tournament size.

    Returns the indexes of the two selected solutions (as a tuple).
    '''
    tour_1 = random.sample(pop_fitness, k=tour_size) # Gets the fitness of a number of random solutions
    sol_index_1 = pop_fitness.index( max(tour_1) ) # Gets the index of the solution with the best fitness in the tournament
    sol_index_2 = sol_index_1
    
    while sol_index_1 == sol_index_2: # This loop condition ensures that the second solution is not the same as the first
        tour_2 = random.sample(pop_fitness, k=tour_size)
        sol_index_2 = pop_fitness.index( max(tour_2) )

    return (sol_index_1, sol_index_2)


def crossover (chrom_1 : list, chrom_2 : list, pc : int) -> list:
    '''
    This function takes two solution chromosomes and either applies the crossover operation on them,
    or copies a random parent, producing an offspring solution.

    :chrom_1: The list representing the first parent's chromosome.
    :chrom_2: The list representing the second parent's chromosome.
    :pc: The crossover rate.

    Returns the offspring's chromosome.
    '''
    offspring = [] # The list for the offspring solution

    if random.uniform(0, 1) <= pc:
        # If the random number between 0 and 1 is less than the crossover rate, then perform crossover.
        for index in range(len(chrom_1)):
            # Selects genes randomly from between the two chromsomes.
            if random.randint(1,2) == 1:
                offspring.append(chrom_1[index])
            else:
                offspring.append(chrom_2[index])
    else:
        # Otherwise, copy the chromosome of one of the parent solutions.
        if random.randint(1,2) == 1:
            offspring = list(chrom_1)
        else:
            offspring = list(chrom_2)
    
    return offspring


def mutate (chrom : list, pm : float, bins : int, which_bpp : bool) -> list:
    '''
    This function takes a solution chromsome and applies the mutation operation to it.

    :chrom: The list representing the chromosome.
    :pm: The mutation rate.
    :bins: The number of bins in the BPP. Defines the domain of the chromosome genes.
    :which_bpp: Set to True for coursework BPP1, or False for BPP2.

    Returns the mutated chromosome.
    '''
    chrom_copy = list(chrom)

    for index in range( len(chrom) ):
        if random.uniform(0, 1) <= pm:
            # If the random number between 0 and 1 is less than the mutation rate, mutate the current gene.
            chrom[index] = random.randint(1, bins)

    if calc_fitness(chrom, bins, which_bpp) < calc_fitness(chrom_copy, bins, which_bpp):
        return chrom_copy
    else:
        return chrom
    
    return chrom


def run (which_bpp : bool, pop_size : int, pm : float, tour_size : int) -> tuple[list, list]:
    '''
    This function runs the genetic algorithm once.

    :which_bpp: Set to True for coursework BPP1, or False for BPP2. Used to determine the number of bins and item weights.
    :pop_size: The size of the population.
    :pm: The mutation rate (probability of mutating each gene in a chromosome).
    :tour_size: The tournament size (number of solutions included in selection tournaments).

    Returns the chromosome of the best solution in the final population, and a list containing the highest fitness in each generation.
    '''
    pc = 0.8 # The crossover rate (probability of crossing over two selected chromosomes)
    elite_size = 1 # The number of best solutions kept from the previous population generation (generational elitism)
    items = 500 # The number of items in this BPP
    bins = 0 # The number of bins in this BPP
    max_fit_evals = 10000 # The maximum number of fitness evaluations before this run of the algorithm is terminated

    if which_bpp == True:
        bins = 10 # 10 bins for BPP1
    else:
        bins = 50 # 50 bins for BPP2

    population = init_population(pop_size, items, bins) # Initialises the population of solutions
    fit_evals = 0 # The number of fitness evaluations
    pop_fitness = [] # The list that will contain the fitness of each solution in the population
    fitness_history = [] # The list that will track the highest fitness value seen in each generation

    while fit_evals < max_fit_evals:
        # Terminate the loop if the maximum number of fitness evaluations has been reached

        pop_fitness = list( map(calc_fitness, population, [bins] * pop_size, [which_bpp] * pop_size) )

        fitness_history.append(max(pop_fitness)) # Adds the highest fitness value in the current population
        
        fit_evals += pop_size # One fitness evaluation for each member of the population, per generation/iteration

        next_gen = [] # Holds the next generation of solutions

        next_gen.append( population[ pop_fitness.index( max(pop_fitness) ) ] )
        # Gets the best solution from the current population and keeps it for the next generation

        while len(next_gen) < pop_size:
            # Create offspring solutions until the correct population size is reached
            
            parents = select_parents(pop_fitness, tour_size) # Indexes of the selected parent solutions
            offspring = crossover(population[parents[0]], population[parents[1]], pc) # Performs crossover on these parent solutions
            offspring = mutate(offspring, pm, bins, which_bpp) # Mutates the offspring solution

            # fit_evals += 2 # The fitness evaluations used in the greedy mutation version

            next_gen.append(offspring)

        population = list(next_gen) # Replaces the current population with the next generation

    pop_fitness = list( map(calc_fitness, population, [bins] * pop_size, [which_bpp] * pop_size) )

    fitness_history.append(max(pop_fitness)) # Adds the final best fitness value to the end of the list.

    best_solution = population[ pop_fitness.index( max(pop_fitness) ) ]
    # Gets the best solution from the final population

    return best_solution, fitness_history


trial_num = 0

while (trial_num < 40):
    # This loop will run 5 trials of the genetic algorithm per combination of mutation rate, tournament size, and BPP (40 total trials)

    random.seed(trial_num) # Seeds the random number generator with a different random number each trial

    pm : float # Mutation rate
    tour_size : int # Tournament size
    which_bpp : bool # True for BPP1, False for BPP2
    bins : int # Number of bins in this BPP

    if (trial_num // 5) % 2 == 0:
        pm = 0.01
    else:
        pm = 0.05

    if (trial_num // 10) % 2 == 0:
        tour_size = 3
    else:
        tour_size = 7

    which_bpp = trial_num // 20 == 0

    if which_bpp == True:
        print(f"BPP1 - Mutation Rate = {pm} - Tournament Size = {tour_size} - Trial {(trial_num % 5) + 1}")
        bins = 10
    else:
        print(f"BPP2 - Mutation Rate = {pm} - Tournament Size = {tour_size} - Trial {(trial_num % 5) + 1}")
        bins = 50

    solution, fitness_history = run(which_bpp, 100, pm, tour_size)
    # Runs the algorithm and returns the chromomsome of the best found solution, along with a list of the best fitness

    bin_weights = [0] * bins

    for gene_index in range(len(solution)):
        # Calculates the weight of each bin in the best solution
        if which_bpp == True: 
            bin_weights[solution[gene_index] - 1] += gene_index + 1 # For BPP1
        else:
            bin_weights[solution[gene_index] - 1] += ((gene_index + 1) ** 2) / 2 # For BPP2

    print("Best Solution Bin Weights : ", bin_weights)
    print(f"Best Solution Fitness : {fitness_history[-1]}")
    # The highest fitness value in the original population vs the highest fitness ever seen.

    try:
        import matplotlib.pyplot as plt # Used to make graphs
        x_axis = range(0, len(fitness_history))
        y_axis = fitness_history
        plt.plot(x_axis, y_axis, linewidth = 8)

        if trial_num % 5 == 4:
            plt.xlabel("Generation", fontsize = 14)  # Label for the X-axis
            plt.ylabel("Fitness", fontsize = 14)  # Label for the Y-axis
            plt.xticks(fontsize=12)
            plt.yticks(fontsize=12)

            if which_bpp == True:
                plt.title(f"BPP1 - Mutation Rate = {pm} - Tournament Size = {tour_size}")  # Chart title
            else:
                plt.title(f"BPP2 - Mutation Rate = {pm} - Tournament Size = {tour_size}")  # Chart title

            plt.show()
    except:
        # If Matplotlib has not been installed
        if (trial_num == 0):
            print("The Matplotlib library has not been installed, so graphs cannot be displayed.")

    print()
    trial_num += 1