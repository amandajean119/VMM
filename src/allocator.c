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

#define MAXGEN 1000
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
#define VCORES 50
#define CROSSRATE 0.2
#define MUTRATE 0.2

/*
 * Check to make sure assignment satisfies formulas
 */

int checkConstraints(pcore * pcore, vm * vm) {
	int tdf = vm->tdf;
	int nrVcores = pcore->nrVCores;
	int uMax = pcore->maxUtilization;
	int pcoreSpeed = pcore->speedKhz;
	int i;
	float utilization = 0.0;
	int speedSum = 0.0;
	float slicePeriod = 0.0;
	for (i = 0; i < nrVcores; i++) {
		slicePeriod = pcore->vcores[i].slice/pcore->vcores[i].period;
		utilization += slicePeriod;
		speedSum += pcore->vcores[i].speedKhz;
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
  int nrVcores = VCORES;
  int i, j, k, satisfies, index;
  pcore * pcores = malloc(nrPcores*sizeof(pcore));
  vm * vm = malloc(sizeof(vm));
  
  for (i = 0; i < nrPcores; i++) {
    srand(time(NULL));
    pcores[i].speedKhz = rand() % config->maxCPU + config->minCPU;
    pcores[i].maxUtilization = rand() % config->maxU + config->minU;
    pcores[i].nrVCores = 0;
    pcores[i].vcores = malloc(nrVcores*sizeof(vcore));
    
  }
  
    vm->vcores = malloc(nrVcores*sizeof(vcore));
    
    for (j = 0;  j < nrVcores; j++) {
      vm->vcores[j].speedKhz =
	rand() % config->maxCPU + config->minCPU;
    }


  for (i = 0; i < ga->populationSize; i++) {
      ga->population[i].pcores = pcores;
      ga->population[i].nrPcores = nrPcores;
      ga->population[i].vm->nrVcores = vm->nrVcores;
      ga->population[i].vm->tdf =
        			rand() % config->maxTdf + config->minTdf;
        for (k = 0; k < nrVcores; k++) {
        	satisfies = -1;
    		ga->population[i].vm->vcores[k].speedKhz =
    				vm->vcores[k].speedKhz;
        	while (satisfies == -1) {
        		ga->population[i].vm->vcores[k].slice =
        				rand() % config->maxSlice + config->minSlice;

        		ga->population[i].vm->vcores[k].period =
						rand() % config->maxPeriod + config->minPeriod;
        		index = rand() % nrPcores;
				ga->population[i].vm->vcores[k].pcore =
					&pcores[index];
				ga->population[i].pcores[index].
					vcores[ga->population[i].pcores[index].nrVCores] =
							ga->population[i].vm->vcores[k];
				satisfies = checkConstraints(ga->population[i].vm->vcores[k].pcore, &ga->population[i].vm);
			}
        	ga->population[i].pcores[index].nrVCores++;
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
 * Perform one-point crossover:
 * Two parents (p1 and p2) are randomly chosen from the current population.
 * Then, a random integer x between 0 and VMS-1 is generated. Two new children
 * (ch1 and ch2) are generated from those parents. The VMs in the positions 0
 * to x of the array of VMs of p1 and the VMs in the positions x+1 to VMS-1 of
 * the array of VMs of p2 are copied to ch1. On the other hand, the rest of VMs
 * of p1 and p2 are copied to ch2. These two new children are added to the new
 * population. This process is repeated until the number of individuals is equal
 * to POPULATION.
 */

/*
 * TODO:
 * CheckConstraints
 */


int crossover(ga * ga){

	int i,j;
	int nrNewPop = 0;
	machine * newPopulation = malloc(POPULATION*sizeof(machine));

	// Two parents chosen randomly

	while(nrNewPop < POPULATION){

		int indexP1 = rand() % ga->populationSize;
		int indexP2 = rand() % ga->populationSize;

		// Parent pointers

		machine * p1 = &ga->population[indexP1];
		machine * p2 = &ga->population[indexP2];

		// Children pointers

		machine * ch1 = malloc(sizeof(machine));
		machine * ch2 = malloc(sizeof(machine));


		// Set the information of the physical cores
		// (It is constant for all machines)

		ch1->nrPcores = PCORES;
		ch1->fitness = p1->fitness; // Should be recalculated
									//if crossover is performed
		ch1->pcores = malloc(PCORES * sizeof(pcore));
		ch1->vm = malloc(sizeof(vm));

		ch2->nrPcores = PCORES;
		ch1->fitness = p2->fitness; // Should be recalculated
									//if crossover is performed
		ch2->pcores = malloc(PCORES * sizeof(pcore));
		ch2->vm = malloc(sizeof(vm));

		for (i = 0; i < PCORES; i++) {

			ch1->pcores[i].maxUtilization = p1->pcores[i].maxUtilization;
			ch1->pcores[i].speedKhz = p1->pcores[i].speedKhz;
			ch1->pcores[i].utilization = 0; // Should be recalculated later

			ch2->pcores[i].maxUtilization = p1->pcores[i].maxUtilization;
			ch2->pcores[i].speedKhz = p1->pcores[i].speedKhz;
			ch2->pcores[i].utilization = 0; // Should be recalculated later

		}
		// Point for crossover generated randomly

		int x = rand() % VCORES;
		if(rand()>CROSSRATE && x > 0){
			ch1->vm->tdf = p1->vm->tdf;
			ch2->vm->tdf = p2->vm->tdf;
			for (i = 0; i < x; i++) {
				ch1->vm->vcores[i].period = p1->vm->vcores[i].period;
				ch1->vm->vcores[i].slice = p1->vm->vcores[i].slice;
				ch1->vm->vcores[i].speedKhz = p1->vm->vcores[i].speedKhz;
				ch1->vm->vcores[i].pcore = malloc(sizeof(pcore));
				ch1->vm->vcores[i].pcore->maxUtilization =
						p1->vm->vcores[i].pcore->maxUtilization;
				ch1->vm->vcores[i].pcore->speedKhz =
						p1->vm->vcores[i].pcore->speedKhz;
				ch1->vm->vcores[i].pcore->utilization =
						p1->vm->vcores[i].pcore->utilization;

				ch2->vm->vcores[i].period = p2->vm->vcores[i].period;
				ch2->vm->vcores[i].slice = p2->vm->vcores[i].slice;
				ch2->vm->vcores[i].speedKhz = p2->vm->vcores[i].speedKhz;
				ch2->vm->vcores[i].pcore = malloc(sizeof(pcore));
				ch2->vm->vcores[i].pcore->maxUtilization =
						p2->vm->vcores[i].pcore->maxUtilization;
				ch2->vm->vcores[i].pcore->speedKhz =
						p2->vm->vcores[i].pcore->speedKhz;
				ch2->vm->vcores[i].pcore->utilization =
						p2->vm->vcores[i].pcore->utilization;
			}
			for (i = x; i < VCORES; i++) {
				ch1->vm->vcores[i].period = p2->vm->vcores[i].period;
				ch1->vm->vcores[i].slice = p2->vm->vcores[i].slice;
				ch1->vm->vcores[i].speedKhz = p2->vm->vcores[i].speedKhz;
				ch1->vm->vcores[i].pcore = malloc(sizeof(pcore));
				ch1->vm->vcores[i].pcore->maxUtilization =
						p2->vm->vcores[i].pcore->maxUtilization;
				ch1->vm->vcores[i].pcore->speedKhz =
						p2->vm->vcores[i].pcore->speedKhz;
				ch1->vm->vcores[i].pcore->utilization =
						p2->vm->vcores[i].pcore->utilization;

				ch2->vm->vcores[i].period = p1->vm->vcores[i].period;
				ch2->vm->vcores[i].slice = p1->vm->vcores[i].slice;
				ch2->vm->vcores[i].speedKhz = p1->vm->vcores[i].speedKhz;
				ch2->vm->vcores[i].pcore = malloc(sizeof(pcore));
				ch2->vm->vcores[i].pcore->maxUtilization =
						p1->vm->vcores[i].pcore->maxUtilization;
				ch2->vm->vcores[i].pcore->speedKhz =
						p1->vm->vcores[i].pcore->speedKhz;
				ch2->vm->vcores[i].pcore->utilization =
						p1->vm->vcores[i].pcore->utilization;
			}
		}
		else{

			ch1->vm->tdf = p1->vm->tdf;
			ch2->vm->tdf = p2->vm->tdf;
			for (j = 0; j < VCORES; ++j) {ch1->vm->vcores[i].period = p1->vm->vcores[i].period;
			ch1->vm->vcores[i].slice = p1->vm->vcores[i].slice;
			ch1->vm->vcores[i].speedKhz = p1->vm->vcores[i].speedKhz;
			ch1->vm->vcores[i].pcore = malloc(sizeof(pcore));
			ch1->vm->vcores[i].pcore->maxUtilization =
					p1->vm->vcores[i].pcore->maxUtilization;
			ch1->vm->vcores[i].pcore->speedKhz =
					p1->vm->vcores[i].pcore->speedKhz;
			ch1->vm->vcores[i].pcore->utilization =
					p1->vm->vcores[i].pcore->utilization;

			ch2->vm->vcores[i].period = p2->vm->vcores[i].period;
			ch2->vm->vcores[i].slice = p2->vm->vcores[i].slice;
			ch2->vm->vcores[i].speedKhz = p2->vm->vcores[i].speedKhz;
			ch2->vm->vcores[i].pcore = malloc(sizeof(pcore));
			ch2->vm->vcores[i].pcore->maxUtilization =
					p2->vm->vcores[i].pcore->maxUtilization;
			ch2->vm->vcores[i].pcore->speedKhz =
					p2->vm->vcores[i].pcore->speedKhz;
			ch2->vm->vcores[i].pcore->utilization =
					p2->vm->vcores[i].pcore->utilization;
			}


		}

		newPopulation[++nrNewPop] = *ch1;
		newPopulation[++nrNewPop] = *ch2;

}

	// TODO CheckConstraints
	free(ga->population);
	ga->population = newPopulation;

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

	int maxGen = MAXGEN;
	int nrGen = 0;			// Number of generations
	struct ga * ga;

	initGA(ga);
	evaluateFitness(ga);


	while(nrGen < maxGen && ga->bestFitness<1.0 ){

		selection(ga);
		crossover(ga);
		mutation(ga);
		evaluateFitness(ga);

	}



}
