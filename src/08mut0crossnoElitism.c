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


#define MAXGEN 50
#define MAXCPU 10000000
#define MINCPU 1000000
#define MAXU 100
#define MINU 80
#define MAXTDF 20
#define MINTDF 1
#define MAXSLICE 1000000
#define MINSLICE 1000
#define MAXPERIOD 5000000
#define MINPERIOD 1000
#define POPULATION 100
#define PCORES 50
#define VCORES 650

#define CROSSRATE 0.9
#define MUTRATE 0.2


/*
 * Free a single machine
 */ 

int freeMachine(machine * m){
 
   int i;	
   printf("Freeing the vcores lists of pcores\n");
   for (i = 0; i < m->nrPcores; i++) {
       printf("Inside for 1\n");
       free(m->pcores[i].vcores);
    }

   printf("Freeing the vcores lists of VM\n");
   for(i=0; i< VCORES; i++){
        printf("Inside for 2\n");
	free(m->vm->vcores[i].pcore);
    }

  printf("Freeing the vcores lists of pcores\n");	
  free(m->pcores);
  printf("Freeing the vcores lists of vcores\n");	
  free(m->vm->vcores);
  printf("Freeing vm\n");	
  free(m->vm);

  return 0;

}

/*
 * Free the whole population
 */


int freePopulation(machine * pop, int popSize){

  //int i;
  //for (i = 0; i < popSize; i++) {
	//freeMachine(&pop[i]);
  //}

  free(pop);
  return 0;

}

/*
 * Free GA Data structure
 */

int freeGA(struct ga * ga) {

  
  freePopulation(ga->population, ga->populationSize);
  free(ga->aConfig);
  free(ga);
  return 0;
}

/*
 * Calculate Utilization of each physical core
 */

float calculateUtilization(pcore * pcore) {

	int nrVcores = pcore->nrVCores;
	int i;
	float utilization = 0.0;
	float slicePeriod = 0.0;
	for (i = 0; i < nrVcores; i++) {
		
		slicePeriod = 100 * (float)pcore->vcores[i].slice/(float)pcore->vcores[i].period;
		utilization += slicePeriod;
	}
	
	return utilization;
}

/*
 * Check to make sure assignment satisfies formulas
 */

int checkConstraints(pcore * pcore, vm * vm) {
	int tdf = vm->tdf;
	int nrVcores = pcore->nrVCores;
	float uMax = pcore->maxUtilization;
	int pcoreSpeed = pcore->speedKhz;
	int i;
	float utilization = calculateUtilization(pcore);
	float speedSum = 0;
	for (i = 0; i < nrVcores; i++) {
		if (pcore->vcores[i].slice > pcore->vcores[i].period) {
			return -1;
		}
	  	speedSum += pcore->vcores[i].speedKhz;
	}
	
	if (utilization/(float)tdf > uMax) {
		return -1;
	}
	if (speedSum/(float)pcoreSpeed > tdf){
		return -1;
	}
	return 0;
}


/*
 * Remove vcore from pcore
 */

int removeVcoreFromPcore(pcore * pc, int vcId){

	if(!pc)	printf("Physical core is NULL\n");	

	else {
		int i;
		int index=-1;
		for(i=0; i < pc->nrVCores; i++){
			if(pc->vcores[i].id == vcId){
				index = i;
				break;			
				}		
		}	
		
                if(index !=-1){
			if(index != (pc->nrVCores-1))
				pc->vcores[index] = pc->vcores[pc->nrVCores-1];

			//vcore * lastVcore = &pc->vcores[pc->nrVCores-1];
			//lastVcore = NULL;	
			pc->nrVCores--;

                }
		else{
			printf("Attempted to remove non existent Vcore\n");
		}// end of else
	} // end of else
	return 0;

}

/*
 * Compute initial population
 */

int initGA(struct ga * ga){

  vcore * vc;
  int index = -1;

 // printf("Initializing general parameters\n");

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
  
  ga->aConfig->maxPeriod = MAXPERIOD;
  ga->aConfig->minPeriod = MINPERIOD;
  ga->aConfig->maxSlice = MAXSLICE;
  ga->aConfig->minSlice = MINSLICE;
 
 
  ga->population = (machine *) malloc(ga->populationSize*sizeof(machine));

  // Configuration which is constant for each machine
  int nrPcores = PCORES;
  int nrVcores = VCORES;
  int i, j, k, satisfies;

  for (i = 0; i < ga->populationSize; i++) {

    ga->population[i].id = i;
    pcore * pcores =  (pcore *) malloc(nrPcores*sizeof(pcore));
    vm * vm = (struct vm *) malloc(sizeof(struct vm));
    
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
      vm->vcores[j].id = j;
    }
    ga->population[i].pcores = pcores;
    ga->population[i].vm = vm;
    ga->population[i].nrPcores = nrPcores;
    ga->population[i].vm->nrVcores = vm->nrVcores;
    
        for (k = 0; k < nrVcores; k++) {
		vc = &ga->population[i].vm->vcores[k];

        	satisfies = -1;
		index = -1;
    		vc->speedKhz = vm->vcores[k].speedKhz;
		vc->id = vm->vcores[k].id;

        	while (satisfies == -1) {
                  
		 if(index!=-1) {
                  removeVcoreFromPcore(&ga->population[i].pcores[index], vc->id);
                 }
		  ga->population[i].vm->tdf =
		    rand() % (ga->aConfig->maxTdf - ga->aConfig->minTdf) + ga->aConfig->minTdf;
		  vc->slice =
		    rand() % (ga->aConfig->maxSlice - ga->aConfig->minSlice) + ga->aConfig->minSlice;
		  
		  vc->period =
		    rand() % (ga->aConfig->maxPeriod - ga->aConfig->minPeriod) + ga->aConfig->minPeriod;
		  index = rand() % nrPcores;
		  vc->pcore =
		    &ga->population[i].pcores[index];
		  ga->population[i].pcores[index].
		    vcores[ga->population[i].pcores[index].nrVCores] = *vc;
		++ga->population[i].pcores[index].nrVCores;
		  satisfies = checkConstraints(vc->pcore, ga->population[i].vm);
		} // end while
		
		ga->population[i].vm->vcores[k] = *vc;
				
		ga->population[i].pcores[index].utilization = calculateUtilization(&(ga->population[i].pcores[index]));
        }
    }
  return 0;
}


/*
 * Evaluate fitnesses of all the individuals of the population
 */

int evaluateFitness(ga * ga){
	float fitness=0.0;
	float utilizationPercent;
	int i, j;
	ga->avgFitness = 0.0;
	ga->bestFitness = 0.0;
	ga->bestMachineIndex = 0;
	for (i = 0; i < ga->populationSize; i++) {
		utilizationPercent = 0.0;
	        fitness = 0.0;
		for (j = 0; j < ga->population[i].nrPcores; j++) {
		  utilizationPercent += ga->population[i].pcores[j].utilization/
		                       (ga->population[i].pcores[j].maxUtilization*ga->population[i].vm->tdf);
		} // end for
		fitness = utilizationPercent/ (float)ga->population[i].nrPcores;
		ga->population[i].fitness = fitness;
		ga->avgFitness += fitness;

		if (fitness >= ga->bestFitness) {
			ga->bestFitness = fitness;
			ga->bestMachineIndex = ga->population[i].id;
		} // end if
	} // end for
	ga->avgFitness /= ga->populationSize;
        return 0;

}

/*
 * Get position in the array population of an individual with identifier = id
 */

int getIndividualIndex(ga * ga, int id){

 int i;
    for(i=0; i<POPULATION;i++){
	if(ga->population[i].id == id){
		 return i;
	}
    }
    printf("Not found\n");
    return -1;

}


/*
 * Get position in the array population of the elite individual
 */


int getEliteIndex(ga * ga){

    int i;
    for(i=0; i<POPULATION;i++){
	if(ga->population[i].id == ga->bestMachineIndex){
		 return i;
	}
    }
    printf("Not found\n");
    return -1;
}


/*
 * Selection function with no elitism
 */

int selectionNoElitism(ga * ga){
  int i, index1, index2;
  machine indiv0 =  ga->population[0];            
  
  machine * newPopulation = (machine*) malloc(POPULATION*sizeof(machine));
  
  for (i = 0; i < POPULATION; i++) {
    do{
      index1 = rand() % POPULATION;
      index2 = rand() % POPULATION;
    }while((index1==0) || (index2==0)); 
    
    if (ga->population[index1].fitness > ga->population[index2].fitness)
      newPopulation[i] = ga->population[index1];
    else
      newPopulation[i] = ga->population[index2];
    //newPopulation[i].id = i;
  }
  
  freePopulation(ga->population, ga->populationSize);
  ga->population = newPopulation;
  return 0;
}

/*
 * Select the best individuals from the population (tournament selection)
 */

int selection(ga * ga){

	int i, index1, index2;
        int id = getEliteIndex (ga);
	machine  elite = ga->population[id];
        machine indiv0 =  ga->population[0];            

	machine * newPopulation = (machine*) malloc(POPULATION*sizeof(machine));
	newPopulation[0] = elite;
        newPopulation[id] = indiv0;

	for (i = 1; i < POPULATION; i++) {
		do{
			index1 = rand() % POPULATION;
			index2 = rand() % POPULATION;
		}while((index1==0) || (index2==0)); 

		if (ga->population[index1].fitness > ga->population[index2].fitness)
			newPopulation[i] = ga->population[index1];
		else
			newPopulation[i] = ga->population[index2];
		//newPopulation[i].id = i;
	}

	freePopulation(ga->population, ga->populationSize);
	ga->population = newPopulation;
	return 0;
}


/*
 * Perform one-point crossover:
 * Two parents (p1 and p2) are randomly chosen from the current population.
 * Then, a random integer x between 0 and VMS-1 is generated. Two new children
 * (ch1 and ch2) are generated from those parents. The VCs in the positions 0
 * to x of the array of VCs of p1 and the VCs in the positions x+1 to VCS-1 of
 * the array of VCs of p2 are copied to ch1. On the other hand, the rest of VCs
 * of p1 and p2 are copied to ch2. These two new children are added to the new
 * population. This process is repeated until the number of individuals is equal
 * to POPULATION.
 */


int crossover(ga * ga){
	
	int i, indexP1, indexP2, pcoreId;
	int crossCounter = 0;
	int nrNewPop = -1;
	int satisfies = -1;
	machine * newPopulation = (machine*) malloc(POPULATION*sizeof(machine));
        newPopulation[0] = ga->population[0]; // elite
	machine * p1;
	machine * p2;
	machine * ch1;
	machine * ch2;

	while(nrNewPop < POPULATION - 1){

		indexP1 = rand() % (ga->populationSize-1) +1;
		indexP2 = rand() % (ga->populationSize-1) +1;

		// Parent pointers

		p1 = &ga->population[indexP1];
		p2 = &ga->population[indexP2];


		//	if((p1->id != newPopulation[0].id) && (p2->id != newPopulation[0].id)){ // avoid elite      
		// Children pointers

		ch1 = (machine *) malloc(sizeof(machine));
		ch2 = (machine *) malloc(sizeof(machine));


		// Set the information of the physical cores
		// (It is constant for all machines)

		ch1->id = nrNewPop+1;		
		ch1->nrPcores = p1->nrPcores;
		ch1->fitness = 0.0;
		ch1->pcores = (pcore *) malloc(ch1->nrPcores * sizeof(pcore));
		ch1->vm = (vm *) malloc(sizeof(vm));
		ch1->vm->vcores = (vcore *) malloc(VCORES*sizeof(vcore));

		ch2->id = nrNewPop+2;
		ch2->nrPcores = p2->nrPcores;
		ch2->fitness = 0.0;
		ch2->pcores = (pcore *) malloc(ch2->nrPcores * sizeof(pcore));
		ch2->vm = (vm *) malloc(sizeof(vm));
		ch2->vm->vcores = (vcore *) malloc(VCORES*sizeof(vcore));
	        
		for (i = 0; i < ch1->nrPcores; i++) {
		    ch1->pcores[i].id = i;
		    ch1->pcores[i].maxUtilization = p1->pcores[i].maxUtilization;
		    ch1->pcores[i].speedKhz = p1->pcores[i].speedKhz;
		    ch1->pcores[i].utilization = 0; // Should be recalculated later
		    ch1->pcores[i].nrVCores = 0;
		    
		    ch2->pcores[i].id = i;
		    ch2->pcores[i].maxUtilization = p2->pcores[i].maxUtilization;
		    ch2->pcores[i].speedKhz = p2->pcores[i].speedKhz;
		    ch2->pcores[i].utilization = 0; // Should be recalculated later
		    ch2->pcores[i].nrVCores = 0;
		

		}

		
		for(i=0;i<VCORES;i++){
			ch1->vm->vcores[i].pcore =  (pcore *) malloc( sizeof(pcore));
			ch2->vm->vcores[i].pcore =  (pcore *) malloc( sizeof(pcore));
		} // end for

		/*for(i=0;i<ch1->nrPcores;i++){
			ch1->pcores[i].vcores = (vcore *) malloc(VCORES*sizeof(vcore));
			ch2->pcores[i].vcores = (vcore *) malloc(VCORES*sizeof(vcore));
				
		} // end for*/


		while (satisfies < 0 && crossCounter <100 ) {

			for(i=0;i<ch1->nrPcores;i++){
				ch1->pcores[i].nrVCores = 0;
				ch2->pcores[i].nrVCores = 0;
				if(ch1->pcores[i].vcores)
					free(ch1->pcores[i].vcores);
				if(ch1->pcores[i].vcores)
					free(ch2->pcores[i].vcores);
				ch1->pcores[i].vcores = (vcore *) malloc(VCORES*sizeof(vcore));
				ch2->pcores[i].vcores = (vcore *) malloc(VCORES*sizeof(vcore));
				
			} // end for

			int x = rand() % VCORES;// Crossover point

			if(rand()>CROSSRATE && x > 0){

				//printf("***** Doing crossover *****\n");

				ch1->vm->tdf = p1->vm->tdf;
				ch2->vm->tdf = p2->vm->tdf;
				for (i = 0; i < x; i++) {
					vcore * vcoreCh1 = &ch1->vm->vcores[i];
					vcore * vcoreCh2 = &ch2->vm->vcores[i];
		
					vcoreCh1->id = p1->vm->vcores[i].id;
					vcoreCh1->period = p1->vm->vcores[i].period;
					vcoreCh1->slice = p1->vm->vcores[i].slice;
					vcoreCh1->speedKhz = p1->vm->vcores[i].speedKhz;
		
					pcoreId =  p1->vm->vcores[i].pcore->id;
					vcoreCh1->pcore = &ch1->pcores[pcoreId];
					vcoreCh1->pcore->vcores[vcoreCh1->pcore->nrVCores] = *vcoreCh1;
					vcoreCh1->pcore->nrVCores++;
					vcoreCh1->pcore->utilization = calculateUtilization(vcoreCh1->pcore);	
	                 	
					vcoreCh2->id = p2->vm->vcores[i].id;
					vcoreCh2->period = p2->vm->vcores[i].period;
					vcoreCh2->slice = p2->vm->vcores[i].slice;
					vcoreCh2->speedKhz = p2->vm->vcores[i].speedKhz;

					pcoreId =  p2->vm->vcores[i].pcore->id;
					vcoreCh2->pcore = &ch2->pcores[pcoreId];
					vcoreCh2->pcore->vcores[vcoreCh2->pcore->nrVCores] = *vcoreCh2;
					vcoreCh2->pcore->nrVCores++;
					vcoreCh2->pcore->utilization = calculateUtilization(vcoreCh2->pcore);	
			
		
				} // end for
				for (i = x; i < VCORES; i++) {


					vcore * vcoreCh1 = &ch1->vm->vcores[i];
					vcore * vcoreCh2 = &ch2->vm->vcores[i];
			
					vcoreCh1->id = p2->vm->vcores[i].id;
					vcoreCh1->period = p2->vm->vcores[i].period;
					vcoreCh1->slice = p2->vm->vcores[i].slice;
					vcoreCh1->speedKhz = p2->vm->vcores[i].speedKhz;

					pcoreId =  p2->vm->vcores[i].pcore->id;
					vcoreCh1->pcore = &ch1->pcores[pcoreId];
					vcoreCh1->pcore->vcores[vcoreCh1->pcore->nrVCores] = *vcoreCh1;
					vcoreCh1->pcore->nrVCores++;
					vcoreCh1->pcore->utilization = calculateUtilization(vcoreCh1->pcore);	

					vcoreCh2->id = p1->vm->vcores[i].id;
					vcoreCh2->period = p1->vm->vcores[i].period;
					vcoreCh2->slice = p1->vm->vcores[i].slice;
					vcoreCh2->speedKhz = p1->vm->vcores[i].speedKhz;

					pcoreId =  p1->vm->vcores[i].pcore->id;
					vcoreCh2->pcore = &ch2->pcores[pcoreId];
					vcoreCh2->pcore->vcores[vcoreCh2->pcore->nrVCores] = *vcoreCh2;
					vcoreCh2->pcore->nrVCores++;
					vcoreCh2->pcore->utilization = calculateUtilization(vcoreCh2->pcore);	

	
				} // end for
			} // end if
			else{
				//printf("***** No crossover just copying to next generation *****\n");
			
				for (i = 0; i < VCORES; i++) {

					vcore * vcoreCh1 = &ch1->vm->vcores[i];
					vcore * vcoreCh2 = &ch2->vm->vcores[i];
		
					vcoreCh1->id = p1->vm->vcores[i].id;
					vcoreCh1->period = p1->vm->vcores[i].period;
					vcoreCh1->slice = p1->vm->vcores[i].slice;
					vcoreCh1->speedKhz = p1->vm->vcores[i].speedKhz;
		
					pcoreId =  p1->vm->vcores[i].pcore->id;
					vcoreCh1->pcore = &ch1->pcores[pcoreId];
					vcoreCh1->pcore->vcores[vcoreCh1->pcore->nrVCores] = *vcoreCh1;
					vcoreCh1->pcore->nrVCores++;
					vcoreCh1->pcore->utilization = calculateUtilization(vcoreCh1->pcore);	
	                 	
					vcoreCh2->id = p2->vm->vcores[i].id;
					vcoreCh2->period = p2->vm->vcores[i].period;
					vcoreCh2->slice = p2->vm->vcores[i].slice;
					vcoreCh2->speedKhz = p2->vm->vcores[i].speedKhz;
			
					pcoreId =  p2->vm->vcores[i].pcore->id;
					vcoreCh2->pcore = &ch2->pcores[pcoreId];
					vcoreCh2->pcore->vcores[vcoreCh2->pcore->nrVCores] = *vcoreCh2;
					vcoreCh2->pcore->nrVCores++;
					vcoreCh2->pcore->utilization = calculateUtilization(vcoreCh2->pcore);	


				} // end for

				} // end else
		
		
			for (i = 0; i < VCORES; i++) {
				satisfies = 0;
			   	satisfies += checkConstraints(ch1->vm->vcores[i].pcore, ch1->vm);
			   	satisfies += checkConstraints(ch2->vm->vcores[i].pcore, ch2->vm);
			} // end for
		crossCounter++;
	
		}  // end while satisfies
		if (crossCounter >= 100) {
			for (i = 0; i < VCORES; i++) {

					vcore * vcoreCh1 = &ch1->vm->vcores[i];
					vcore * vcoreCh2 = &ch2->vm->vcores[i];
		
					vcoreCh1->id = p1->vm->vcores[i].id;
					vcoreCh1->period = p1->vm->vcores[i].period;
					vcoreCh1->slice = p1->vm->vcores[i].slice;
					vcoreCh1->speedKhz = p1->vm->vcores[i].speedKhz;
		
					pcoreId =  p1->vm->vcores[i].pcore->id;
					vcoreCh1->pcore = &ch1->pcores[pcoreId];
					vcoreCh1->pcore->vcores[vcoreCh1->pcore->nrVCores] = *vcoreCh1;
					vcoreCh1->pcore->nrVCores++;
					vcoreCh1->pcore->utilization = calculateUtilization(vcoreCh1->pcore);	
	                 	
					vcoreCh2->id = p2->vm->vcores[i].id;
					vcoreCh2->period = p2->vm->vcores[i].period;
					vcoreCh2->slice = p2->vm->vcores[i].slice;
					vcoreCh2->speedKhz = p2->vm->vcores[i].speedKhz;

					pcoreId =  p2->vm->vcores[i].pcore->id;
					vcoreCh2->pcore = &ch2->pcores[pcoreId];
					vcoreCh2->pcore->vcores[vcoreCh2->pcore->nrVCores] = *vcoreCh2;
					vcoreCh2->pcore->nrVCores++;
					vcoreCh2->pcore->utilization = calculateUtilization(vcoreCh2->pcore);	


				}
		}

		newPopulation[++nrNewPop] = *ch1;
                if(nrNewPop < POPULATION-1)
			newPopulation[++nrNewPop] = *ch2;
		satisfies = -1;
		crossCounter = 0;
		//	    } // end if elite
        } // end while population
        
        freePopulation(ga->population, ga->populationSize);
	ga->population = newPopulation;
	return 0;

} // end crossover function


/*
 * Mutate individuals
 */

int mutation(ga * ga){

	int i, j, mutAlloc, oldPeriod, oldSlice, newAllocationIndex;
	machine * indiv;
	vcore * vc;
	int oldPcoreIndex;
	pcore * oldPcore;
	pcore * newPcore;

	int satisfies, counter;
	for (i = 1; i < POPULATION; i++) {
	    indiv = &ga->population[i];

	    //	    if(indiv->id != ga->bestMachineIndex){ // Do not mutate elite 
        
	        for (j = 0; j < VCORES; j++) {
		      satisfies = -1;
		      counter = 0;
		      mutAlloc = -1;
		      while (satisfies == -1 && counter < 100) {
		        
		        mutAlloc = -1;
			vc = &indiv->vm->vcores[j];
			oldPeriod = vc->period;
			oldSlice = vc->slice;
		
			if(rand() > MUTRATE){
			  vc->period =
			    rand() % (ga->aConfig->maxPeriod - ga->aConfig->minPeriod) + ga->aConfig->minPeriod;
			
			}
			if(rand() > MUTRATE){
			  vc->slice =
			    rand() % (ga->aConfig->maxSlice - ga->aConfig->minSlice) + ga->aConfig->minSlice;
			}
		
			oldPcoreIndex = vc->pcore->id;
			oldPcore = &indiv->pcores[oldPcoreIndex];
		
			if(rand() > MUTRATE){
			  
			  // We must add vc to the list of vcores of vc->pcore
			  // and remove it from the list of vcores of the old pcore and decrement the nrVcores in that pcore
			  mutAlloc = 0;
			  removeVcoreFromPcore(oldPcore, vc->id);   
			  newAllocationIndex = (rand() % indiv->nrPcores);
			  vc->pcore = &indiv->pcores[newAllocationIndex];
			  
			  vc->pcore->vcores[vc->pcore->nrVCores] = *vc;
			  indiv->pcores[newAllocationIndex].nrVCores++;
			  vc->pcore->utilization = calculateUtilization(vc->pcore);
		          newPcore = &indiv->pcores[newAllocationIndex];
			  
			} //end if
			satisfies = checkConstraints(vc->pcore, ga->population[i].vm);
			if (satisfies == -1){
			  vc->period = oldPeriod;
			  vc->slice = oldSlice;
			  if (mutAlloc == 0) {
		 	    vc->pcore = oldPcore;	
			    oldPcore->vcores[oldPcore->nrVCores] = *vc;
			    oldPcore->nrVCores++;
			    newPcore->nrVCores--;
			    newPcore->utilization = calculateUtilization(newPcore);
			    oldPcore->utilization = calculateUtilization(oldPcore);
			  } // end if
		
			} // end if
			counter++;
		      } //end while
		    } //end for vcores
		//     } // end of if elite
	} // end for population
	
	return 0;
}


/*
 * Print chromosome parameter values
 */

int printChromosome(machine * m, FILE * chromosomeFile){

	int i,j;
	for (i = 0; i < m->nrPcores; i++) {   
		for (j = 0; j <  m->pcores[i].nrVCores ; ++j) { 
			fprintf(chromosomeFile,"%d, %f, %d, %d, %f, %f, %d, %f, %f, %d\n", 
				m->id, 
				m->fitness,
				m->vm->tdf,
				m->pcores[i].id,
				m->pcores[i].utilization/m->vm->tdf, 
				m->pcores[i].maxUtilization, 
				m->pcores[i].vcores[j].id,
				(float)m->pcores[i].vcores[j].slice,
				(float)m->pcores[i].vcores[j].period,
				m->pcores[i].vcores[j].pcore->id);
		} //end for
		
	} // end for

	return 0;
}

/*
 * Print parameters of a machine (used for debugging purposes)
 */


int printMachineParameters(machine * m){

	int i,j;
	printf("********* Machine Parameters *********\n");
	printf("Machine number: %d, Fitness: %f, TDF: %d, #Pcores %d \n", 
		m->id, 
		m->fitness,
		m->vm->tdf,
		m->nrPcores);

	printf("MACHINE %d PHYSICAL CORES:\n", m->id);
	for (i = 0; i < m->nrPcores; i++) {
		printf("PHYSICAL CORE: %d, MaxU: %f, U:%f, #Vcores: %d\n", 
			m->pcores[i].id, 
			m->pcores[i].maxUtilization, 
			m->pcores[i].utilization, 
			m->pcores[i].nrVCores);

		printf("PHYSICAL CORE %d VIRTUAL CORES:\n", m->pcores[i].id);

		for (j = 0; j <  m->pcores[i].nrVCores ; ++j) {
			printf("Vcore: %d, Slice: %f, Period: %f, CPUSpeed: %d, Pcore %d\n",
				m->pcores[i].vcores[j].id,
				(float)m->pcores[i].vcores[j].slice,
				(float)m->pcores[i].vcores[j].period,
				m->pcores[i].vcores[j].speedKhz,
				m->pcores[i].vcores[j].pcore->id);
		}

	}

	return 0;
}

/*
 * Main function
 */

int main(){

	//printf(" *** Real time multiprocessor allocator***\n");

	int i;
	int maxGen = MAXGEN+1;
	int nrGen = 0;			// Number of generations
	ga * ga = (struct ga *) malloc(sizeof(struct ga));
	
	FILE * fitnessFile;
	FILE * chromosomeFile;

	remove("fitness08mut0crossnoElitism.txt");
	remove("chromosome08mut0crossnoElitism.txt");
	
	fitnessFile = fopen("fitness08mut0crossnoElitism.txt", "w");
	chromosomeFile = fopen("chromosome08mut0crossnoElitism.txt", "w");
	
	fprintf(fitnessFile, "%s", "Number of generations, Best Machine, Best fitness, Average fitness \n");
	fprintf(chromosomeFile, "%s", "Machine number, Fitness, TDF, PHYSICAL CORE id, U, Umax, Vcore id, Slice, Period, Pcore\n");
	
	srand(time(NULL));

	printf("Initializing GA\n");
	initGA(ga);
	printf("Evaluating fitness\n");
	evaluateFitness(ga);
        
	fprintf(fitnessFile, "%d, %d, %f, %f\n", nrGen, ga->bestMachineIndex,ga->bestFitness, ga->avgFitness);
	fprintf(chromosomeFile, "%s", "Generation 0\n");
	for (i = 0; i < ga->populationSize; i++){ 
		    if(ga->population[i].id == ga->bestMachineIndex)	
	            	printChromosome(&(ga->population[i]), chromosomeFile);
		}
	
	nrGen++;
	
	while(nrGen < maxGen && ga->bestFitness<1.0 ){

         	printf("Selection for gen %d ...\n", nrGen);
		//selection(ga);
		selectionNoElitism(ga);
		printf("Mutation for gen %d ...\n", nrGen);
		mutation(ga);
        	printf("Crossover for gen %d ...\n", nrGen);
		crossover(ga);
		printf("Evaluating fitness for gen %d ...\n", nrGen);
        	evaluateFitness(ga);
	
        	fprintf(fitnessFile, "%d, %d, %f, %f\n", nrGen, ga->bestMachineIndex, ga->bestFitness, ga->avgFitness);
		fprintf(chromosomeFile, "Generation %d\n", nrGen);
		printChromosome(&(ga->population[0]), chromosomeFile);
        
		nrGen++;
			
	} // end while
 
	fclose(chromosomeFile);
	fclose(fitnessFile);
	
        freeGA(ga);
	return 0;

}
