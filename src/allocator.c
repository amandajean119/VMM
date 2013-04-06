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

#define MAXCPU 10000000
#define MINCPU 1000000
#define MAXU 100
#define MINU 80
#define MAXTDF 20
#define MINTDF 1
#define MAXSLICE 1000000
#define MINSLICE 1000
#define MAXPERIOD 100000000
#define MINPERIOD 100000
#define POPULATION 100
#define PCORES 10
#define VMS 1
#define VCORES 50

/*
 * Check to make sure assignment satisfies formulas
 */

int checkConstraints(pcore * pcore, vm * vm) {
	int tdf = vm->tdf;
	int nrVcores = vm->nrVcores;
	int uMax = pcore->maxUtilization;
	int pcoreSpeed = pcore->speedKhz;
	int i;
	float utilization = 0.0;
	int speedSum = 0.0;
	float slicePeriod = 0.0;
	for (i = 0; i < nrVcores; i++) {
		slicePeriod = vm->vcores[i].slice/vm->vcores[i].period;
		utilization += slicePeriod;
		speedSum += vm->vcores[i].speedKhz;
	}
	if (utilization/tdf > uMax)
		return -1;
	if (speedSum/pcoreSpeed > tdf)
		return -1;
	return 0;
}

/*
 * Compute initial population
 */

int initGA(struct ga * ga){

  allocConfig * config = malloc(sizeof(allocConfig));
  
  // CPU Speeds in Khz
  config->maxCPU = MAXCPU;   // Khz
  config->minCPU = MINCPU ;
  
  // Percentage of utilization
  config->maxU = MAXU;
  config->minU = MINU;
  
  // Time dilation Factor
  config->maxTdf = MAXTDF;
  config->minTdf = MINTDF;
  
  // Slices and periods in micro seconds
  config->maxSlice = MAXSLICE; // Usec
  config->minSlice = MINSLICE;
  config->maxPeriod = MAXPERIOD;
  config->minPeriod = MINPERIOD;
  
  ga->aConfig = config;
  ga->populationSize = POPULATION;
  ga->population = malloc(ga->populationSize*sizeof(machine));
  
  // Configuration which is constant for each machine
  int nrPcores = PCORES;
  int nrVms = VMS;
  int nrVcores = VCORES;
  int i, j, k, l, satisfies;
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
        	satisfies = -1;
    		ga->population[i].vms[j].vcores[k].speedKhz =
    				vms[j].vcores[k].speedKhz;
        	while (satisfies == -1) {
        		ga->population[i].vms[j].vcores[k].slice =
        				rand() % config->maxSlice + config->minSlice;

        		ga->population[i].vms[j].vcores[k].period =
						rand() % config->maxPeriod + config->minPeriod;

				ga->population[i].vms[j].vcores[k].pcore =
					&pcores[rand() % nrPcores];

				satisfies = checkConstraints(ga->population[i].vms[j].vcores[k].pcore, &ga->population[i].vms[j]);
			}
        }
      }
    }

  return 0;
}


/*
 * Evaluate fitnesses of all the individuals of the population
 */

int evaluateFitness(ga * ga){
	float fitness;
	float utilizationPercent;
	int i, j;
	ga->avgFitness = 0.0;
	ga->bestFitness = 0.0;
	for (i = 0; i < ga->populationSize; ++i) {
		utilizationPercent = 0;
		for (j = 0; j < ga->population[i].nrPcores; ++j) {
			utilizationPercent += ga->population[i].pcores[j].utilization/
					ga->population[i].pcores[j].maxUtilization;
		}
		fitness = utilizationPercent/ ga->population[i].nrPcores;
		ga->population[i].fitness = fitness;
		ga->avgFitness += fitness;
		if (fitness > ga->bestFitness)
			ga->bestFitness = fitness;
	}
	ga->avgFitness /= ga->populationSize;
	return 0;

}

/*
 * Select the best individuals from the population
 */

int selection(ga * ga){

	return 0;

}


/*
 * Perform crossover with selected individuals
 */

int crossover(ga * ga){

	return 0;

}

/*
 * Mutate individuals
 */

int mutation(ga * ga){

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
