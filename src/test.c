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


#define MAXGEN 100
#define MAXCPU 10000000
#define MINCPU 1000000
#define MAXU 100
#define MINU 80
#define MAXTDF 20
#define MINTDF 1
#define MAXSLICE 1000000
#define MINSLICE 1000
#define MAXPERIOD 4000000
#define MINPERIOD 1000
#define POPULATION 100
//#define PCORES 50
//#define VCORES 650
#define PCORES 5
#define VCORES 10

#define CROSSRATE 0.2
#define MUTRATE 1.0


int freeMachine(machine * m){
 /*
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
*/
  return 0;

}

int freePopulation(machine * pop, int popSize){

  int i;
  for (i = 0; i < popSize; i++) {
	freeMachine(&pop[i]);
  }

  free(pop);
  return 0;

}

int freeGA(struct ga * ga) {

  
  freePopulation(ga->population, ga->populationSize);
  free(ga->aConfig);
  free(ga);
  return 0;
}


float calculateUtilization(pcore * pcore) {

	int nrVcores = pcore->nrVCores;
	int i;
	float utilization = 0.0;
	float slicePeriod = 0.0;
	for (i = 0; i < nrVcores; i++) {
		
//	printf("Calculate utilization. vcore %d Slice: %f Period: %f\n",pcore->vcores[i].id, 
	//	(float)pcore->vcores[i].slice, 
	//	(float)pcore->vcores[i].period);

		
		slicePeriod = 100 * (float)pcore->vcores[i].slice/(float)pcore->vcores[i].period;
		utilization += slicePeriod;
	}
	//printf("Max Utilization: %f, Actual Utilization from calcUt function: %f\n", pcore->maxUtilization, utilization);
	
	return utilization;
}

/*
 * Check to make sure assignment satisfies formulas
 */

int checkConstraints(pcore * pcore, vm * vm) {
	int tdf = vm->tdf;
	int nrVcores = pcore->nrVCores;
	float uMax = pcore->maxUtilization;
	//printf("in checkConstraints max utilization for pcore %d: %f\n", pcore->id, pcore->maxUtilization);
	int pcoreSpeed = pcore->speedKhz;
	int i;
	float utilization = calculateUtilization(pcore);
	float speedSum = 0;
	//printf("Number of virtual cores for pcore %d: %d\n", pcore->id, nrVcores);
	for (i = 0; i < nrVcores; i++) {
		if (pcore->vcores[i].slice > pcore->vcores[i].period) {
			//printf("****************REJECTED DUE TO SLICE > PERIOD********\n");
			return -1;
		}
	  	speedSum += pcore->vcores[i].speedKhz;
	}
	
	
    	/*printf("********** check constrains ***********\n");
	 printf("check constrains, Pcore id %d\n",pcore->id);
	 printf("check constrains, tdf %d\n",tdf);
	 printf("check constrains, uMax %f\n",uMax);
	 printf("check constrains, speedSum %f\n",speedSum);
	 printf("check constrains, utilization %f\n",utilization);
	 printf("check constrains, pcoreSpeed %d\n",pcoreSpeed);
         printf("check constrains, utilization/tdf  %f\n",utilization/tdf );
	 printf("check constrains, speedSum/pcoreSpeed  %f\n", speedSum/pcoreSpeed );
	 printf("***************************************\n");
      	*/ 
	if (utilization/(float)tdf > uMax) {
		//printf("*********REJECTED DUE TO UTILIZATION*****\n");
		return -1;
	}
	if (speedSum/(float)pcoreSpeed > tdf){
		//printf("*********REJECTED DUE TO SPEED SUM*****\n");
		return -1;
	}
	//printf("************ACCEPTED********\n");
	return 0;
}

int getnrVCores(pcore * pcore) {
  int nrVCores = 0;
  while(&(pcore->vcores[nrVCores]) != NULL)
    ++nrVCores;
  return nrVCores;
      
}

/*
 * Remove vcore from pcore
 */

int removeVcoreFromPcore(pcore * pc, int vcId){

//	printf("***** removeVcoreFromPcore *****\n");
	if(!pc)	printf("********** pc is NULL\n");	

	else {
		//printf("********** pc id %d\n", pc->id);		
		//printf("********** vc id %d\n", vcId);			
		int i;
		int index=-1;
		for(i=0; i < pc->nrVCores; i++){
			//printf("for %d\n",i);			
			//printf("pc->vcores[i].id %d\n",pc->vcores[i].id);
			//printf("pc->vcores[i].period %f\n",(float)pc->vcores[i].period);
			//printf("vcId %d\n",vcId);

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

  //printf("Initializing before for\n");

  for (i = 0; i < ga->populationSize; i++) {

 //   
    ga->population[i].id = i;
//printf("next machine %d\n", ga->population[i].id);
    pcore * pcores =  (pcore *) malloc(nrPcores*sizeof(pcore));
    vm * vm = (struct vm *) malloc(sizeof(struct vm));
    
    //printf("Initializing physical cores\n");
    
    for (j = 0; j < nrPcores; j++) {
      pcores[j].id = j;
      pcores[j].vcores = (vcore *) malloc(nrVcores*sizeof(vcore));
      pcores[j].nrVCores = 0;
      
      pcores[j].maxUtilization = rand() % (ga->aConfig->maxU - ga->aConfig->minU) + ga->aConfig->minU;
   //   printf("max utilization for pcore %d: %f\n", j, pcores[j].maxUtilization);
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
//		printf("for machine %d\n", ga->population[i].id);
//		printf("max utilization for pcore %d: %f\n", index, ga->population[i].pcores[index].maxUtilization);
//printf("max utilization for pcore %d: %f\n", index, vc->pcore->maxUtilization);
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
	float fitness;
	float utilizationPercent;
	int i, j;
	ga->avgFitness = 0.0;
	//ga->bestFitness = 0.0;
	ga->bestMachineIndex = 0;
	for (i = 0; i < ga->populationSize; ++i) {
	//	printf("for machine %d\n", i);
		utilizationPercent = 0.0;
		for (j = 0; j < ga->population[i].nrPcores; ++j) {
		  utilizationPercent += ga->population[i].pcores[j].utilization/
		    (ga->population[i].pcores[j].maxUtilization*ga->population[i].vm->tdf);
	//	printf("utilization percent for pcore %d: %f\n", j, ga->population[i].pcores[j].utilization/
		//    (ga->population[i].pcores[j].maxUtilization*ga->population[i].vm->tdf));
		}
		fitness = utilizationPercent/ (float)ga->population[i].nrPcores;
	//	printf("fitness: %f\n", fitness);
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

int selectionNoElitism(ga * ga){
  //	machine elite = ga->population[ga->bestMachineIndex];
	machine * newPopulation = (machine*) malloc(POPULATION*sizeof(machine));
	//	newPopulation[0] = elite;
	//newPopulation[0].id = 0;
	ga->bestMachineIndex = 0;
	int i, index1, index2;
	for (i = 0; i < POPULATION; i++) {
		index1 = rand() % POPULATION;
		index2 = rand() % POPULATION;
		//printf("Fitness1 %f\n", ga->population[index1].fitness);
		//printf("Fitness2 %f\n", ga->population[index2].fitness);
		if (ga->population[index1].fitness > ga->population[index2].fitness)
			newPopulation[i] = ga->population[index1];
		else
			newPopulation[i] = ga->population[index2];
		newPopulation[i].id = i;
	}
	free(ga->population);
	ga->population = newPopulation;
	return 0;
}

/*
 * Select the best individuals from the population
 */

int selection(ga * ga){
	machine elite = ga->population[ga->bestMachineIndex];
	machine * newPopulation = (machine*) malloc(POPULATION*sizeof(machine));
	newPopulation[0] = elite;
	newPopulation[0].id = 0;
	ga->bestMachineIndex = 0;
	int i, index1, index2;

	for (i = 1; i < POPULATION; i++) {
		index1 = rand() % POPULATION;
		index2 = rand() % POPULATION;
		//printf("Fitness1 %f\n", ga->population[index1].fitness);
		//printf("Fitness2 %f\n", ga->population[index2].fitness);
		if (ga->population[index1].fitness > ga->population[index2].fitness)
			newPopulation[i] = ga->population[index1];
		else
			newPopulation[i] = ga->population[index2];
		newPopulation[i].id = i;
	}

	freePopulation(ga->population, ga->populationSize);
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
	
	int i, indexP1, indexP2, pcoreId;
	int crossCounter = 0;
	int nrNewPop = -1;
	int satisfies = -1;
	machine * newPopulation = (machine*) malloc(POPULATION*sizeof(machine));
	machine * p1;
	machine * p2;
	machine * ch1;
	machine * ch2;

	while(nrNewPop < POPULATION - 1){

		indexP1 = rand() % ga->populationSize;
		indexP2 = rand() % ga->populationSize;

		// Parent pointers

		p1 = &ga->population[indexP1];
		p2 = &ga->population[indexP2];
               
		// Children pointers

	        ch1 = (machine *) malloc(sizeof(machine));
	        ch2 = (machine *) malloc(sizeof(machine));


		// Set the information of the physical cores
		// (It is constant for all machines)

		ch1->id = p1->id;		
		ch1->nrPcores = p1->nrPcores;
		ch1->fitness = 0.0;
		ch1->pcores = (pcore *) malloc(ch1->nrPcores * sizeof(pcore));
		ch1->vm = (vm *) malloc(sizeof(vm));
		ch1->vm->vcores = (vcore *) malloc(VCORES*sizeof(vcore));

		ch2->id = p2->id;
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

			//printf("Inside crossover while\n");
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
		newPopulation[++nrNewPop] = *ch2;
		satisfies = -1;
		crossCounter = 0;
 		
   
	} // end while population

        freePopulation(ga->population, ga->populationSize);
	ga->population = newPopulation;
	return 0;

} // end crossover function


/*
 * Mutate individuals
 */

int mutation(ga * ga){
	

  	int i, j, mutAlloc, oldPeriod, newPeriod, oldSlice, newSlice, newAllocationIndex;
	machine * indiv;
	vcore * vc;
	int oldPcoreIndex;
	pcore * oldPcore;
	pcore * newPcore;

	/*for (i = 0; i < POPULATION; i++) {
		for(j =0; j< VCORES; j++){
			printf("machine %d, vcore %d\n",i,ga->population[i].vm->vcores[j].id);
		} // end for

	}*/ // end for

	int satisfies, counter;
	//printf("IN MUTATION\n");
	for (i = 0; i < POPULATION; ++i) {
	  // satisfies = -1;
	 
	  indiv = &ga->population[i];
	 
	    /*	  oldTDF = indiv->vm->tdf;
		  if(rand() > MUTRATE){
		  newTdf = rand() % (ga->aConfig->maxTdf - ga->aConfig->minTdf) + ga->aConfig->minTdf;
		  //newTdf = indiv->vm->tdf + (rand()%2 -1);
		  if(newTdf>= MINTDF && newTdf<= MAXTDF)
		  indiv->vm->tdf = newTdf;
		  }
	    */
	    for (j = 0; j < VCORES; j++) {
	      satisfies = -1;
	      counter = 0;
	      mutAlloc = -1;
	      while (satisfies == -1 && counter < 100) {
                
                mutAlloc = -1;
		//	printf("inside mutation while, individual %d, vcore %d\n", i, j);
	     
		vc = &indiv->vm->vcores[j];
		oldPeriod = vc->period;
		oldSlice = vc->slice;
		
		//printf("IN MUTATION WHILE LOOP\n");
		if(rand() > MUTRATE){
		  vc->period =
		    rand() % (ga->aConfig->maxPeriod - ga->aConfig->minPeriod) + ga->aConfig->minPeriod;
		  //newPeriod = vc->period + (rand() % 2 - 1);
		  //if(newPeriod>= MINPERIOD && newPeriod<= MAXPERIOD) 
		  //vc->period = newPeriod;

		}
		if(rand() > MUTRATE){
		  vc->slice =
		    rand() % (ga->aConfig->maxSlice - ga->aConfig->minSlice) + ga->aConfig->minSlice;
		  //		  newSlice = vc->slice + (rand() % 2 - 1);
		  //if(newSlice>= MINSLICE && newSlice<= MAXSLICE)
		  //vc->slice = newSlice;
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
	} // end for population
	
	return 0;
}


int printChromosome(machine * m, FILE * chromosomeFile){

	int i,j;
	for (i = 0; i < m->nrPcores; i++) {
		for (j = 0; j <  m->pcores[i].nrVCores ; ++j) {
			fprintf(chromosomeFile,"%d, %d, %d, %f, %d, %f, %f, %d\n", 
				m->id, 
				m->vm->tdf,
				m->pcores[i].id,
				m->pcores[i].utilization, 
				m->pcores[i].vcores[j].id,
				(float)m->pcores[i].vcores[j].slice,
				(float)m->pcores[i].vcores[j].period,
				m->pcores[i].vcores[j].pcore->id);
		} //end for
		
	} // end for

	return 0;
}


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
	
	remove("test_fitness.txt.txt");
	remove("test_chromosome.txt");
	
	fitnessFile = fopen("test_fitness.txt", "w");
	chromosomeFile = fopen("test_chromosome.txt", "w");
	
	fprintf(fitnessFile, "%s", "Number of generations, Best Machine, Best fitness, Average fitness \n");
	fprintf(chromosomeFile, "%s", "Machine number, TDF, PHYSICAL CORE id, U, Vcore id, Slice, Period, Pcore\n");
	
	srand(time(NULL));

	//printf("Initializing GA\n");
	initGA(ga);
	//printf("Evaluating fitness\n");
	evaluateFitness(ga);
	fprintf(fitnessFile, "%d, %d, %f, %f\n", nrGen, ga->bestMachineIndex,ga->bestFitness, ga->avgFitness);
	fprintf(chromosomeFile, "%s", "Generation 0\n");
	for (i = 0; i < ga->populationSize; i++) 
	  printChromosome(&(ga->population[i]), chromosomeFile);
	
	
	nrGen++;
	//for (i = 0; i < ga->populationSize; i++) 
	//	printMachineParameters(&(ga->population[i]));

	//while(nrGen < maxGen && ga->bestFitness<1.0 ){
	while(nrGen < maxGen && ga->bestFitness<1.0 ){

		printf("Selection for gen %d ...\n", nrGen);
		selection(ga);
		printf("Mutation for gen %d ...\n", nrGen);
		mutation(ga);
		printf("Crossover for gen %d ...\n", nrGen);
		crossover(ga);
		printf("Evaluating fitness for gen %d ...\n", nrGen);
		evaluateFitness(ga);
		fprintf(fitnessFile, "%d, %d, %f, %f\n", nrGen, ga->bestMachineIndex, ga->bestFitness, ga->avgFitness);
		fprintf(chromosomeFile, "Generation %d", nrGen);		
		for (i = 0; i < ga->populationSize; i++) 
		  printChromosome(&(ga->population[i]), chromosomeFile);
	
	/*	printf("Number of generations %d, ", nrGen);
		printf("Best fitness: %f, ", ga->bestFitness);
		printf("Average fitness: %f \n", ga->avgFitness);

		//for (i = 0; i < ga->populationSize; i++)
		//  printMachineParameters(&(ga->population[i]));
	*/
	
		nrGen++;
			
	} // end while
 
	fclose(chromosomeFile);
	fclose(fitnessFile);
	
        freeGA(ga);
	return 0;

}
