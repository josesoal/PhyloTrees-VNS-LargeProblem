/*
 ============================================================================
 Authors :
 	Jose Luis Soncco-Alvarez and Mauricio Ayala-Rincon
 	Group of Theory of Computation
 	Universidade de Brasilia (UnB) - Brazil
 ============================================================================
 */

#ifndef DCJ_DIST_H_
#define DCJ_DIST_H_

int DCJdistance( PointDCJPtr *genome1DCJ, PointDCJPtr *genome2DCJ, 
				int *inverseGenome1, int *inverseGenome2, 
				int numPoints1DCJ, int numPoints2DCJ, int numberGenes );
void calculateInverseGenome( PointDCJPtr *genomeDCJ, int numPointsDCJ, int *inverseGenome );
void applyDCJ( PointDCJPtr *genomeDCJ, int *numPointsDCJ, int i, int j, int firstForm );


#endif

