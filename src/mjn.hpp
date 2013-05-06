#ifndef _MJN_H
#define _MKN_H

#include <Rcpp.h>
#include <algorithm> 
#include <iostream> 
#include <igraph/igraph.h>

using namespace Rcpp;

typedef std::pair<int, int> Combination;

// some debugging print functions
void printDataVectors (std::vector<std::vector<int> > x);
void printCombos (std::set<Combination> combos);
void printIntVector (std::vector<int> x);

// helper functions
std::set<Combination> allPairs (int n);
int nNodes (std::vector<std::vector<int> > *edges);

// R to C++ helper
std::vector<std::vector<int> > convert2Matrix (IntegerMatrix *x);
std::vector<int> convert2Vector (IntegerVector *x);

// core functionality
std::vector<std::vector<int> > manDist(std::vector<std::vector<int> > *x);
std::vector<std::vector<int> > hammingDist(std::vector<std::vector<int> > *x, std::vector<int> *w);

std::vector<bool> connected (std::vector<std::vector<int> > *edges, 
                             igraph_t *g, igraph_vector_t *weights, 
                             int delta, int epsilon);

void addFeasibleEdges (std::vector<std::vector<int> > *edges, 
                       std::vector<int> *deltas, int epsilon, int delta,
                       igraph_t *g, igraph_vector_t *weights, 
                       std::vector<bool> *feasible);

std::vector<std::vector<int> > medianVectors(std::vector<std::vector<int> > triplet);

void eraseDuplicates (std::vector<std::vector<int> > *new_vectors, 
                      std::vector<std::vector<int> > *existing);

std::vector<int> calcConnectionCosts(std::vector<std::vector<int> > triplet, 
                                     std::vector<std::vector<int> > new_vectors,
                                     std::vector<int> *character_weights);

int newSequenceTypes (std::vector<std::vector<int> > *data,
                  std::vector<int> *character_weights,
                  std::vector<std::vector<int> > *edges,
                  std::vector<bool> feasible,
                  int epsilon,
                  std::vector<std::vector<int> > *mvecs);
                  
std::vector<bool> deltaStepComponents (std::vector<std::vector<int> > *edges,
                                       std::vector<int> *deltas,
                                       std::set<int> *unique_deltas,
                                       int epsilon);

std::set<int, std::greater<int> > identifyObsoletes (std::vector<std::vector<int> > *edges,
                                                     std::vector<bool> *feasible,
                                                     int num_sampled,
                                                     int num_total);

void mergeAndConvertData(std::vector<std::vector<int> > *alldata, 
                         std::vector<std::vector<int> > *data, 
                         std::vector<std::vector<int> > *mvecs, 
                         std::vector<int> *deltas, 
                         std::set<int> *unique_deltas,
                         std::vector<std::vector<int> > *edges);
                         
void mjnC (std::vector<std::vector<int> > *data,
           std::vector<int> *character_weights,
           std::vector<std::vector<int> > *edgeList, 
           std::vector<int> *edgeLengths, 
           int epsilon);
           
// interface exported to R
List mjnRcpp (IntegerMatrix dataR, IntegerVector charWeights, int epsilon);

#endif