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
#define POPULATION 2
#define PCORES 10
#define VCORES 30
#define CROSSRATE 0.2
#define MUTRATE 0.2


float calculateUtilization(pcore * pcore) {

	int nrVcores = pcore->nrVCores;
	printf("nrVCores: %d\n", nrVcores);
	int i;
	float utilization = 0.0;
	float slicePeriod = 0.0;
	for (i = 0; i < nrVcores; i++) {
		slicePeriod = 100 * (float)pcore->vcores[i].slice/(float)pcore->vcores[i].period;
		utilization += slicePeriod;
		printf("current utilization %f\n", utilization);
	}
	printf("done with this pcore\n");
	return utilization;
}

/*
 * Check to make sure assignment satisfies formulas
 */

int checkConstraints(pcore * pcore, vm * vm) {
	int tdf = vm->tdf;
	int nrVcores = pcore->nrVCores;
	int uMax = pcore->maxUtilization;
	int pcoreSpeed = pcore->speedKhz;
	int i;
	float utilization = calculateUtilization(pcore);
	int speedSum = 0;
	float slicePeriod = 0.0;
	for (i = 0; i < nrVcores; i++) {
	  speedSum += pcore->vcores[i].speedKhz;
	}

	//	printf("check constrains, tdf %d\n",tdf);
	//	printf("check constrains, uMax %d\n",uMax);
	//  printf("check constrains, speedSum %d\n",speedSum);
	// printf("check constrains, pcoreSpeed %d\n",pcoreSpeed);

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

  //printf("Initializing general parameters\n");

   ga-> aConfig  =  (struct allocatorConfig *) malloc(sizeof(struct allocatorConfig));
 
  // CPU Speeds in Khz
  ga->aConfig->maxCPU = MAXCPU;   // Khz
  ga->aConfig->minCPU = MINCPU ;

  // Percentage of utilization
  ga->aConfig->maxU = MAXU;
  ga->aConfig->minU = MINU;

  // Time dilation Factor
  ga->aConfig->maxTdf = MAXTDF;
  ga->aConfig->minTdf = MINTDF;
  ga->populationSize = POPULATION;
  // Slices and periods in micro seconds
  //WEIRD BUG, SEG FAULTS IF PERIOD PARAMS ARE BELOW SLICE PARAMS
  ga->aConfig->maxPeriod = MAXPERIOD;
  ga->aConfig->minPeriod = MINPERIOD;
  ga->aConfig->maxSlice = MAXSLICE; // Usec
  ga->aConfig->minSlice = MINSLICE;
 
 
  ga->population = (machine *) malloc(ga->populationSize*sizeof(machine));

  // Configuration which is constant for each machine
  int nrPcores = PCORES;
  int nrVcores = VCORES;
  int i, j, k, satisfies, index;

  //printf("Initializing before for\n");

  for (i = 0; i < ga->populationSize; i++) {

    //    printf("next machine %d\n", i);

	  pcore * pcores =  (pcore *) malloc(nrPcores*sizeof(pcore));
	  vm * vm = (struct vm *) malloc(sizeof(struct vm));

	  //	  printf("Initializing physical cores\n");

	  for (j = 0; j < nrPcores; j++) {
	    pcores[j].id = j;
	    pcores[j].vcores = (vcore *) malloc(nrVcores*sizeof(vcore));
	    pcores[j].nrVCores = 0;

	    pcores[j].maxUtilization = rand() % (ga->aConfig->maxU - ga->aConfig->minU) + ga->aConfig->minU;

	    pcores[j].speedKhz = rand() % (ga->aConfig->maxCPU - ga->aConfig->minCPU) + ga->aConfig->minCPU;
	    
	  }

	    vm->vcores = (struct vcore *) malloc(nrVcores*sizeof(struct vcore));

	    for (j = 0;  j < nrVcores; j++) {
	      vm->vcores[j].speedKhz =
		(rand() % (ga->aConfig->maxCPU - ga->aConfig->minCPU)) + ga->aConfig->minCPU;
	      //printf("vcore speed: %d\n", vm->vcores[j].speedKhz);
	      vm->vcores[j].id = j;
	    }
      ga->population[i].pcores = pcores;
      ga->population[i].vm = vm;
      ga->population[i].nrPcores = nrPcores;
      ga->population[i].vm->nrVcores = vm->nrVcores;

        for (k = 0; k < nrVcores; k++) {
        	satisfies = -1;
    		ga->population[i].vm->vcores[k].speedKhz =
    				vm->vcores[k].speedKhz;
        	while (satisfies == -1) {

		  //	printf("Inside while\n");
        	      ga->population[i].vm->tdf =
        	        			rand() % (ga->aConfig->maxTdf - ga->aConfig->minTdf) + ga->aConfig->minTdf;
        		ga->population[i].vm->vcores[k].slice =
        				rand() % (ga->aConfig->maxSlice - ga->aConfig->minSlice) + ga->aConfig->minSlice;

        		ga->population[i].vm->vcores[k].period =
						rand() % (ga->aConfig->maxPeriod - ga->aConfig->minPeriod) + ga->aConfig->minPeriod;
        		index = rand() % nrPcores;
				ga->population[i].vm->vcores[k].pcore =
					&pcores[index];
				ga->population[i].pcores[index].
					vcores[ga->population[i].pcores[index].nrVCores] =
							ga->population[i].vm->vcores[k];
				satisfies = checkConstraints(ga->population[i].vm->vcores[k].pcore, ga->population[i].vm);
			}
		++ga->population[i].pcores[index].nrVCores;
		ga->population[i].pcores[index].utilization = calculateUtilization(&(ga->population[i].pcores[index]));
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
	ga->bestMachineIndex = 0;
	for (i = 0; i < ga->populationSize; ++i) {
		utilizationPercent = 0;
		for (j = 0; j < ga->population[i].nrPcores; ++j) {
			utilizationPercent += ga->population[i].pcores[j].utilization/
					ga->population[i].pcores[j].maxUtilization;
		}
		fitness = utilizationPercent/ ga->population[i].nrPcores;
		ga->population[i].fitness = fitness;
		ga->avgFitness += fitness;
		if (fitness > ga->bestFitness) {
			ga->bestFitness = fitness;
			ga->bestMachineIndex = i;
		}
	}
	ga->avgFitness /= ga->populationSize;
	return 0;

}

/*
 * Select the best individuals from the population
 */

int selection(ga * ga){
	machine elite = ga->population[ga->bestMachineIndex];
	machine * newPopulation = malloc(POPULATION*sizeof(machine));
	newPopulation[0] = elite;
	ga->bestMachineIndex = 0;
	int i, index1, index2;
	for (i = 1; i < POPULATION; i++) {
		index1 = rand() % POPULATION;
		index2 = rand() % POPULATION;
		if (ga->population[index1].fitness > ga->population[index2].fitness)
			newPopulation[i] = ga->population[index1];
		else
			newPopulation[i] = ga->population[index2];
	}
	free(ga->population);
	ga->population = newPopulation;
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


int crossover(ga * ga){

  int i, j, indexP1, indexP2;
	int nrNewPop = 0;
	int satisfies = -1;
	machine * newPopulation = malloc(POPULATION*sizeof(machine));
	machine * p1;
	machine * p2;
	machine * ch1;
	machine * ch2;
	// Two parents chosen randomly

	while(nrNewPop < POPULATION){

		indexP1 = rand() % ga->populationSize;
		indexP2 = rand() % ga->populationSize;

		// Parent pointers

		p1 = &ga->population[indexP1];
		p2 = &ga->population[indexP2];

		// Children pointers

	        ch1 = malloc(sizeof(machine));
	        ch2 = malloc(sizeof(machine));


		// Set the information of the physical cores
		// (It is constant for all machines)

		ch1->nrPcores = PCORES;
		ch1->fitness = p1->fitness; // Should be recalculated
									//if crossover is performed
		ch1->pcores = malloc(ch1->nrPcores * sizeof(pcore));
		ch1->vm = malloc(sizeof(vm));

		ch2->nrPcores = PCORES;
		ch2->fitness = p2->fitness; // Should be recalculated
									//if crossover is performed
		ch2->pcores = malloc(ch2->nrPcores * sizeof(pcore));
		ch2->vm = malloc(sizeof(vm));
		for (i = 0; i < ch1->nrPcores; i++) {

			ch1->pcores[i].maxUtilization = p1->pcores[i].maxUtilization;
			ch1->pcores[i].speedKhz = p1->pcores[i].speedKhz;
			ch1->pcores[i].utilization = 0; // Should be recalculated later

			ch2->pcores[i].maxUtilization = p2->pcores[i].maxUtilization;
			ch2->pcores[i].speedKhz = p2->pcores[i].speedKhz;
			ch2->pcores[i].utilization = 0; // Should be recalculated later

		}
		// Point for crossover generated randomly
		while (satisfies < 0) {
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
				for (j = 0; j < VCORES; ++j) {
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


			}
			satisfies = 0;
			for (i = 0; i < VCORES; i++) {
				satisfies += checkConstraints(ch1->vm->vcores[i].pcore, ch1->vm);
				satisfies += checkConstraints(ch2->vm->vcores[i].pcore, ch2->vm);
			}
		}

		newPopulation[++nrNewPop] = *ch1;
		newPopulation[++nrNewPop] = *ch2;
		for (i = 0; i < newPopulation[nrNewPop].nrPcores; i++) {
		  //WHAT ABOUT NUMBER OF VCORES??
		  newPopulation[nrNewPop-1].pcores[i].utilization = calculateUtilization(&(newPopulation[nrNewPop-1].pcores[i]));
		  newPopulation[nrNewPop].pcores[i].utilization = calculateUtilization(&(newPopulation[nrNewPop].pcores[i]));
										      
		}
		
	}

	free(ga->population);
	ga->population = newPopulation;

	return 0;
}

/*
 * Mutate individuals
 */

int mutation(ga * ga){

	int i, j;
	int satisfies = -1;

	for (i = 0; i < POPULATION; ++i) {
		machine * indiv = &ga->population[i];

		if(rand() > MUTRATE){
			int newTdf = indiv->vm->tdf + (rand()%2 -1);
			if(newTdf>= MINTDF && newTdf<= MAXTDF)
				indiv->vm->tdf = newTdf;
		}

		for (j = 0; j < VCORES; ++j) {
			while (satisfies == -1) {
				if(rand() > MUTRATE){
					int newPeriod = indiv->vm->vcores[j].period + (rand() % 2 - 1);
					if(newPeriod>= MINPERIOD && newPeriod<= MAXPERIOD)
						indiv->vm->vcores[i].period = newPeriod;
				}
				if(rand() > MUTRATE){
					int newSlice = indiv->vm->vcores[j].slice + (rand() % 2 - 1);
					if(newSlice>= MINSLICE && newSlice<= MAXSLICE)
						indiv->vm->vcores[i].slice = newSlice;
				}
				if(rand() > MUTRATE){
					int newAllocationIndex = (rand() % indiv->nrPcores);
					indiv->vm->vcores[i].pcore = &indiv->pcores[newAllocationIndex];
				}
			satisfies = checkConstraints(ga->population[i].vm->vcores[j].pcore, ga->population[i].vm);
			}
			++ga->population[i].pcores[j].nrVCores;
			ga->population[i].pcores[j].utilization = calculateUtilization(&(ga->population[i].pcores[j]));
		}
	}

	return 0;
}

int printMachineParameters(machine * m){

	int i,j;
	//	pcore pCore;
	//	vcore vCore;
	printf("Fitness %f\n", m->fitness);
	printf("Time dilation factor %d\n",	m->vm->tdf);
	printf("Number of pcores %d\n", m->nrPcores);
	
	for (i = 0; i < m->nrPcores; i++) {

	  //pCore = m->pcores[i];
		printf("For pcore number %d:\n", m->pcores[i].id);
		printf("Maximum utilization %d\n",  m->pcores[i].maxUtilization);
		printf("Actual utilization %d\n",  m->pcores[i].utilization);
		printf("Number of virtual cores %d \n",  m->pcores[i].nrVCores);
		printf("For the virtual cores allocated to pcore number %d: \n", m->pcores[i].id);


		for (j = 0; j <  m->pcores[i].nrVCores ; ++j) {

		  //vCore =  m->pcores[i].vcores[j];
			printf("For vcore number %d:\n",m->pcores[i].vcores[j].id);
		    printf("Slice %f\n",(float)m->pcores[i].vcores[j].slice);
		    printf("Period %f",(float)m->pcores[i].vcores[j].period);
		    printf("CPU Speed (Khz) %d \n",m->pcores[i].vcores[j].speedKhz);
            printf("Physical core id %d \n", m->pcores[i].vcores[j].pcore->id);

		}

	}

	return 0;
}
/*
int free(ga) {
  int i, j, k, l;
  for (i = 0; i < ga->populationSize; i++) {
    for (j = 0; j < ga->machine[i]->nrPcores; j++) {
      for (k = 0; k < ga->machine[i]->pcores[j].nrVCores; k++) {
	free(ga->machine[i]->pcores[j].vcore[k]);
      }
    }
    for (j = 0; j < ga-> machine[i]->vm->nrVcores; j++) {
      free(ga->machine[i]->vm->vcores[j]);
  }
  }
  return 0;
}
*/
/*
 * Main function
 */

int main(){

	printf(" *** Real time multiprocessor allocator***\n");

	int i;
	int maxGen = MAXGEN;
	int nrGen = 0;			// Number of generations
	ga * ga = (struct ga *) malloc(sizeof(struct ga));
	srand(time(NULL));

	
	initGA(ga);

	/*	
	evaluateFitness(ga);

	for (i = 0; i < ga->populationSize; i++) {
		printMachineParameters(&(ga->population[i]));

	}

	


	while(nrGen < maxGen && ga->bestFitness<1.0 ){

		selection(ga);
		crossover(ga);
		mutation(ga);
		evaluateFitness(ga);

	}*/

        free(ga);
	return 0;

}
