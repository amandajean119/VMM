/*
 * allocator.c
 *
 *  Created on: Apr 3, 2013
 *  Author:
 */

#include "allocator.h"

/*
 * Compute initial population
 */

int initGA(struct ga * ga){

	return 0;

}

/*
 * Evaluate fitnesses of all the individuals of the population
 */

int evaluateFitness(struct ga * ga){

	return 0;

}

/*
 * Select the best individuals from the population
 */

int selection(struct ga * ga){

	return 0;

}


/*
 * Perform crossover with selected individuals
 */

int crossover(struct ga * ga){

	return 0;

}

/*
 * Mutate individuals
 */

int mutation(struct ga * ga){

	return 0;

}


int maxFitness(struct ga * ga){

	return 0;

}



/*
 * Main function
 */

int main(){

	printf(" *** Real time multiprocessor allocator***\n");

	int maxGen = 1000;
	int nrGen = 0;			// Number of generations
	struct ga * ga;

	initGA(ga);
	evaluateFitness(ga);

	/*
	while(nrGen < maxGen && maxFitness(ga)<100 ){

		selection(ga);
		crossover(ga);
		mutation(ga);
		evaluateFitness(ga);

	}
 */


}
