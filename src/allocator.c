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
#define VMS 1
#define VCORES 50
#define CROSSRATE 0.2
#define MUTRATE 0.2

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

		// Parents pointers

		machine * p1 = &ga->population[indexP1];
		machine * p2 = &ga->population[indexP2];

		// Children pointers

		machine * ch1 = malloc(sizeof(machine));
		machine * ch2 = malloc(sizeof(machine));


		// Set the information of the physical cores
		// (It is constant for all machines)

		ch1->nrPcores = PCORES;
		ch1->nrVms = VMS;
		ch1->fitness = p1->fitness; // Should be recalculated
									//if crossover is performed
		ch1->pcores = malloc(PCORES * sizeof(pcore));
		ch1->vms = malloc(VMS * sizeof(vm));

		ch2->nrPcores = PCORES;
		ch2->nrVms = VMS;
		ch1->fitness = p2->fitness; // Should be recalculated
									//if crossover is performed
		ch2->pcores = malloc(PCORES * sizeof(pcore));
		ch2->vms = malloc(VMS * sizeof(vm));

		for (i = 0; i < PCORES; i++) {

			ch1->pcores[i].maxUtilization = p1->pcores[i].maxUtilization;
			ch1->pcores[i].speedKhz = p1->pcores[i].speedKhz;
			ch1->pcores[i].utilization = 0; // Should be recalculated later

			ch2->pcores[i].maxUtilization = p1->pcores[i].maxUtilization;
			ch2->pcores[i].speedKhz = p1->pcores[i].speedKhz;
			ch2->pcores[i].utilization = 0; // Should be recalculated later

		}

		if(rand()>CROSSRATE){

			// Point for crossover generated randomly

			int x = rand() % VMS;


			// According with the crossover point (x) the VMS of the parents are
			// assigned to the children.

			for (i = 0; i < VMS; i++) {

				vm * vmCh1 = &ch1->vms[i];
				vmCh1->nrVcores = VCORES;
				vmCh1->vcores = malloc(VCORES*sizeof(vcore));

				vm * vmCh2 = &ch2->vms[i];
				vmCh2->nrVcores = VCORES;
				vmCh2->vcores = malloc(VCORES*sizeof(vcore));


				if(i <= x){

					vmCh1->nrVcores = p1->vms[i].nrVcores;
					vmCh1->tdf = p1->vms[i].tdf;

					for (j = 0; j < VCORES; j++) {

						vmCh1->vcores[j].pcore = p1->vms[i].vcores[j].pcore;
						vmCh1->vcores[j].period = p1->vms[i].vcores[j].period;
						vmCh1->vcores[j].slice = p1->vms[i].vcores[j].slice;
						vmCh1->vcores[j].speedKhz = p1->vms[i].vcores[j].speedKhz;

					}

					vmCh2->nrVcores = p2->vms[i].nrVcores;
					vmCh2->tdf = p2->vms[i].tdf;

					for (j = 0; j < VCORES; j++) {

						vmCh2->vcores[j].pcore = p2->vms[i].vcores[j].pcore;
						vmCh2->vcores[j].period = p2->vms[i].vcores[j].period;
						vmCh2->vcores[j].slice = p2->vms[i].vcores[j].slice;
						vmCh2->vcores[j].speedKhz = p2->vms[i].vcores[j].speedKhz;

					}
				}

				else {

					vmCh1->nrVcores = p2->vms[i].nrVcores;
					vmCh1->tdf = p2->vms[i].tdf;

					for (j = 0; j < VCORES; j++) {

						vmCh1->vcores[j].pcore = p2->vms[i].vcores[j].pcore;
						vmCh1->vcores[j].period = p2->vms[i].vcores[j].period;
						vmCh1->vcores[j].slice = p2->vms[i].vcores[j].slice;
						vmCh1->vcores[j].speedKhz = p2->vms[i].vcores[j].speedKhz;

					}

					vmCh2->nrVcores = p1->vms[i].nrVcores;
					vmCh2->tdf = p1->vms[i].tdf;

					for (j = 0; j < VCORES; j++) {

						vmCh2->vcores[j].pcore = p1->vms[i].vcores[j].pcore;
						vmCh2->vcores[j].period = p1->vms[i].vcores[j].period;
						vmCh2->vcores[j].slice = p1->vms[i].vcores[j].slice;
						vmCh2->vcores[j].speedKhz = p1->vms[i].vcores[j].speedKhz;

					}

				}

			}
		}

		else{

			for (i = 0; i < VMS; ++i) {

				ch1->vms[i].nrVcores = p1->vms[i].nrVcores;
				ch1->vms[i].tdf = p1->vms[i].tdf;

				ch2->vms[i].nrVcores = p2->vms[i].nrVcores;
				ch2->vms[i].tdf = p2->vms[i].tdf;

				for (j = 0; j < VCORES; ++j) {

					ch1->vms[i].vcores[j].period = p1->vms[i].vcores[j].period;
					ch1->vms[i].vcores[j].slice = p1->vms[i].vcores[j].slice;
					ch1->vms[i].vcores[j].speedKhz = p1->vms[i].vcores[j].speedKhz;
					ch1->vms[i].vcores[j].pcore = malloc(sizeof(pcore));
					ch1->vms[i].vcores[j].pcore->maxUtilization =
							ch1->vms[i].vcores[j].pcore->maxUtilization;
					ch1->vms[i].vcores[j].pcore->speedKhz =
							ch1->vms[i].vcores[j].pcore->speedKhz;
					ch1->vms[i].vcores[j].pcore->utilization =
							ch1->vms[i].vcores[j].pcore->utilization;

					ch2->vms[i].vcores[j].period = p2->vms[i].vcores[j].period;
					ch2->vms[i].vcores[j].slice = p2->vms[i].vcores[j].slice;
					ch2->vms[i].vcores[j].speedKhz = p2->vms[i].vcores[j].speedKhz;
					ch2->vms[i].vcores[j].pcore = malloc(sizeof(pcore));
					ch2->vms[i].vcores[j].pcore->maxUtilization =
							ch2->vms[i].vcores[j].pcore->maxUtilization;
					ch2->vms[i].vcores[j].pcore->speedKhz =
							ch2->vms[i].vcores[j].pcore->speedKhz;
					ch2->vms[i].vcores[j].pcore->utilization =
							ch2->vms[i].vcores[j].pcore->utilization;


				}

			}



		}

		newPopulation[++nrNewPop] = *ch1;
		newPopulation[++nrNewPop] = *ch2;

	}


	for (i = 0; i < PCORES; ++i) {

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
