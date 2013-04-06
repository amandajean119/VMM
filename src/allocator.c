/*
 * allocator.c
 *
 *  Created on: Apr 3, 2013
 *  Author:
 */

#include "allocator.h"
#include <stdio.h>      /* printf, scanf, NULL */
#include <stdlib.h>     /* malloc, free, rand */
#include <time.h>       /* srand */

/*
 * Compute initial population
 */

int initGA(struct ga * ga){

  allocConfig * config = malloc(sizeof(allocConfig));
  
  // CPU Speeds in Khz
  config->maxCPU = 12000000;   // Khz
  config->minCPU = 300000;
  
  // Percentage of utilization
  config->maxU = 100;
  config->minU = 80;
  
  // Time dilation Factor
  config->maxTdf = 20;
  config->minTdf = 1;
  
  // Slices and periods in micro seconds
  config->maxSlice = 1000000; // Usec
  config->minSlice = 1000;
  config->maxPeriod = 100000000;
  config->minPeriod = 100000;
  
  ga->aConfig = config;
  ga->populationSize = 100;
  ga->population = malloc(ga->populationSize*sizeof(machine));
  
  // Configuration which is constant for each machine
  int nrPcores = 10;
  int nrVms = 1;
  int nrVcores = 50;
  int i, j, k, l;
  pcore * pcores = malloc(nrPcores*sizeof(pcore));
  vm * vms = malloc(nrVms*sizeof(vm));
  
  for (i = 0; i < nrPcores; i++) {
    srand(time(NULL));
    pcores[i].speedKhz = rand() % config->maxCPU + config->minCPU;
    pcores[i].maxUtilization = rand() % config->maxU + config->minU;
    
  }
  
  for (i = 0; i < nrVms; i++) {
    
    vms[i].vcores = malloc(nrVcores*sizeof(vcore));
    
    for (j = 0;  j < nrVcores; j++) {
      vms[i].vcores[j].speedKhz =
	rand() % config->maxCPU + config->minCPU;
    }
  }

  for (i = 0; i < ga->populationSize; i++) {
      ga->population[i].pcores = pcores;
      ga->population[i].nrPcores = nrPcores;
      ga->population[i].nrVms = nrVms;

      for (j = 0;  j < nrVms; j++) {
        ga->population[i].vms[j].nrVcores = vms[j].nrVcores;
        ga->population[i].vms[j].tdf =
        			rand() % config->maxTdf + config->minTdf;

        for (k = 0; k < nrVcores; k++) {

        	ga->population[i].vms[j].vcores[k].slice =
					rand() % config->maxSlice + config->minSlice;

			ga->population[i].vms[j].vcores[k].period =
								rand() % config->maxPeriod + config->minPeriod;

			ga->population[i].vms[j].vcores[k].speedKhz =
					vms[j].vcores[k].speedKhz;
			while (!equationissatisfied) {
				ga->population[i].vms[j].vcores[k].pcore =
					&pcores[rand() % nrPcores];
			}
        }
      }
    }

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


	while(nrGen < maxGen && maxFitness(ga)<100 ){

		selection(ga);
		crossover(ga);
		mutation(ga);
		evaluateFitness(ga);

	}



}
