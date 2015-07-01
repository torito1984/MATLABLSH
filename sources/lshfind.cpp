// Adaptation of LSH for euclidean spaces for MATLAB
// Author: David Martinez Rego

#include <stdio.h>
#include <stdlib.h>

#include "mex.h"
#include "headers.h"

#define N_SAMPLE_QUERY_POINTS 100

// Total memory free in GB (Required by the original library)
// Set high for avoiding its definition 15 GB (can end up wihtout memory)
// If large scale, set this value higher
// Other possibility is to let it go until it runs out of memory
// By changing the definition of MALLOC

#define TOTAL_AVAILABLE_MEMORY 16106127360

void populatePoints(int dimension, int nData, double* data, PPointT* dataPoints)
{
    int i, j;
    for(i = 0; i < nData; i++) {
        FAILIF(NULL == (dataPoints[i] = (PPointT)MALLOC(sizeof(PointT))));
        dataPoints[i]->index = i;
        FAILIF(NULL == (dataPoints[i]->coordinates = (RealT*)MALLOC(dimension * sizeof(RealT))));
        for(j = 0; j < dimension; j++) {
            dataPoints[i]->coordinates[j] = data[(i * dimension) + j];
        }
    }
}

void generateRandomQueries(int nSampleQueries, int nData, PPointT* dataPoints, PPointT* sampleQueries)
{
    for(IntT i = 0; i < nSampleQueries; i++) {
        sampleQueries[i] = dataPoints[genRandomInt(0, nData - 2)];
    }
}

void cleanPoints(int nData, PPointT* dataPoints)
{
    for(int i = 0; i < nData; i++) {
        FREE(dataPoints[i]->coordinates)
        FREE(dataPoints[i])
    }
    FREE(dataPoints)
}

void cleanup(int nData, PPointT* dataPoints, int nQueries, PPointT* queryPoints, PPointT* sampleQueries){
    cleanPoints(nData, dataPoints);
    cleanPoints(nQueries, queryPoints);
    FREE(sampleQueries);
}


void searchPoints(double radius, double succP, size_t dimensions, size_t nQueries, size_t nData, 
                double *queries, double *data, double maxReported, int nlhs, mxArray* plhs[]){
     
    // For result management
    int nnFound;
    IntT resultSize;
    PPointAndRealTStructT *distToNN = NULL;
    int cPos = 0;
    double* out;
    int *tmp;
   
    // LSH data struct
    PPointT* dataPoints = NULL;
    PPointT* queryPoints = NULL;
    PPointT* sampleQueries = NULL;
    PRNearNeighborStructT nnStruct = NULL;
    RNNParametersT optParameters;
    
     // Construct the vector of points
    FAILIF(NULL == (dataPoints = (PPointT*)MALLOC(nData * sizeof(PPointT))));
    populatePoints(dimensions, nData, data, dataPoints);
    FAILIF(NULL == (queryPoints = (PPointT*)MALLOC(nQueries * sizeof(PPointT))));
    populatePoints(dimensions, nQueries, queries, queryPoints);

    // Calculate the optimal parameters in memory
    FAILIF(NULL == (sampleQueries = (PPointT*)MALLOC(N_SAMPLE_QUERY_POINTS * sizeof(PPointT))));
    generateRandomQueries(N_SAMPLE_QUERY_POINTS, nData, dataPoints, sampleQueries);
    
    optParameters = computeOptimalParameters(radius,
                                             succP,
                                             nData,
                                             dimensions,
                                             dataPoints,
                                             N_SAMPLE_QUERY_POINTS,
                                             sampleQueries,
                                             (MemVarT)((availableTotalMemory - totalAllocatedMemory)));
                                           
        
    // Calculate neighbors
    nnStruct = initLSH_WithDataSet(optParameters, nData, dataPoints);
 
    resultSize = nData;
    PPointT* result = (PPointT*)MALLOC(resultSize * sizeof(PPointT));

    tmp = (int *)MALLOC((maxReported+2)*nQueries*sizeof(int));
    cPos = 0;
   
    for(IntT i = 0; i < nQueries; i++) {
        nnFound = getRNearNeighbors(nnStruct, queryPoints[i], result, resultSize);
        if(nnFound > 0) {
            // compute the distances to the found NN, and sort according to the distance
            FAILIF(NULL == (distToNN = (PPointAndRealTStructT*)REALLOC(distToNN, nnFound * sizeof(PPointAndRealTStructT))));
            for(IntT p = 0; p < nnFound; p++) {
                distToNN[p].ppoint = result[p];
                distToNN[p].real = distance(dimensions,  queryPoints[i], result[p]);
            }
            qsort(distToNN, nnFound, sizeof(PPointAndRealTStructT), comparePPointAndRealTStructT);

            tmp[cPos] = i+1; cPos++; // query
            tmp[cPos] = MIN(nnFound, maxReported); cPos++; // numFound
            
            // Out the NNs 
            for(IntT j = 0; j < MIN(nnFound, maxReported); j++) {
               tmp[cPos] = distToNN[j].ppoint->index+1; cPos++;
            }
            FREE(distToNN);
        }
    }
    
    // allocate and initialise output matrix
    plhs[0] = mxCreateDoubleMatrix(1, cPos, mxREAL);
    out  = mxGetPr(plhs[0]);
    for(IntT j = 0; j < cPos; j++)
        out[j] = tmp[j];
            
   
    // Clear. we do not want leaks!!
    FREE(result); 
    FREE(tmp);
    freePRNearNeighborStruct(nnStruct); 
    cleanup(nData, dataPoints, nQueries, queryPoints, sampleQueries);
}

// [res] = lshfind( R, succP, queries, set, maxReported);
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]){

    double radius, succP, maxReported;

    // Pointers to the real data
    double* queries;
    double* data;
    
    size_t dimensions;
    size_t nQueries;
    size_t nData;

    radius = mxGetPr(prhs[0])[0];
    succP = mxGetPr(prhs[1])[0];

    queries = mxGetPr(prhs[2]);
    data = mxGetPr(prhs[3]);

    dimensions = mxGetM(prhs[2]);
    nQueries = mxGetN(prhs[2]);
    nData = mxGetN(prhs[3]);
    maxReported = mxGetPr(prhs[4])[0];

    // Register available memory
    totalAllocatedMemory = 0;
    availableTotalMemory = TOTAL_AVAILABLE_MEMORY;
    
    searchPoints(radius, succP, dimensions, nQueries, nData, queries, data, maxReported, nlhs, plhs);

}

/**************************** That's all Folks!  *****************************/
