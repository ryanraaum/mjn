#include "mjn.hpp"
// #include <ctime>

// time_t starttime, endtime;

void 
printDataVectors (std::vector<std::vector<int> > x)
{
  for (std::vector<std::vector<int> >::iterator xit=x.begin(); xit!=x.end(); ++xit) {
    for (std::vector<int>::iterator yit=xit->begin(); yit!=xit->end(); ++yit) {
      std::cout << *yit << " ";
    }
    std::cout << std::endl;
  }  
}

void 
printCombos (std::set<Combination> combos)
{
  std::set<Combination>::iterator cit;
  
  for (cit=combos.begin(); cit!=combos.end(); ++cit) {
    std::cout << cit->first << " " << cit->second << std::endl;
  }
}


void 
printIntVector (std::vector<int> x)
{
  for (std::vector<int>::iterator it=x.begin(); it!=x.end(); ++it) {
    std::cout << *it << " ";
  }  
  std::cout << std::endl;
}


std::vector<std::vector<int> > 
manDist(std::vector<std::vector<int> > *x) 
{
    int sum;
    int nrows = x->size();
    int ncols = (*x)[0].size();
    std::vector<std::vector<int> > D (nrows, std::vector<int> (nrows));
    
    for (int i = 0; i < (nrows-1); i++) {
      for (int j = i + 1; j < nrows; j++) {
        sum = 0;
        for (int k = 0; k < ncols; k++) {
          sum = sum + abs((*x)[i][k] - (*x)[j][k]);
        }
        D[i][j] = sum;
        D[j][i] = sum;
      }
    }
    return D;
}

std::vector<std::vector<int> > 
hammingDist(std::vector<std::vector<int> > *x, std::vector<int> *w) 
{
    int sum;
    int nrows = x->size();
    int ncols = (*x)[0].size();
    std::vector<std::vector<int> > D (nrows, std::vector<int> (nrows));
    
    for (int i = 0; i < (nrows-1); i++) {
      for (int j = i + 1; j < nrows; j++) {
        sum = 0;
        for (int k = 0; k < ncols; k++) {
          if ((*x)[i][k] != (*x)[j][k]) {
            sum = sum + (*w)[k];
          }
        }
        D[i][j] = sum;
        D[j][i] = sum;
      }
    }
    return D;
}


void 
connected (std::vector<std::vector<int> > *edges, 
           std::vector<int> *deltas,
           igraph_t *g, 
           igraph_vector_t *weights, 
           int delta, 
           int epsilon,
           std::vector<bool> *already_connected) 
{
  // std::vector<bool> result (edges->size(), false);
  igraph_bool_t areconnected;
  igraph_vector_t epath, vpath;
  int nedges;
  int eindex; 
  long int numedges;
  bool all_less;
  
  //std::cout << "number of edges: " << edges->size() << std::endl;
  for (int i = 0; i < edges->size(); i++) {
    if ((((*deltas)[i] - epsilon) <= delta) && !(*already_connected)[i]) {
      igraph_st_edge_connectivity(g, &nedges, (*edges)[i][0], (*edges)[i][1]);
      igraph_are_connected(g, (*edges)[i][0], (*edges)[i][1], &areconnected);
      
      if (nedges == 0) {
        (*already_connected)[i] = false;
      } else if ((nedges > 0) && areconnected) {
        (*already_connected)[i] = true;
      } else if (nedges > 0) {
        igraph_vector_init(&vpath, 0);
        igraph_vector_init(&epath, 0);
        igraph_get_shortest_path_dijkstra(g, &vpath, &epath, (*edges)[i][0], (*edges)[i][1], weights, IGRAPH_ALL);
        numedges = igraph_vector_size(&epath);
        all_less = true;
        for (int k = 0; k < numedges; k++) {
          eindex = VECTOR(epath)[k];
          all_less = all_less && (VECTOR(*weights)[eindex] < (delta - epsilon));
        }
        (*already_connected)[i] = all_less;
        igraph_vector_destroy(&epath);
        igraph_vector_destroy(&vpath);
      }
    }
  }
  // return result;
}

int 
nNodes (std::vector<std::vector<int> > *edges) 
{
  std::vector<int> all_nodes;
  for (int i=0; i < edges->size(); i++) {
    for (int j=0; j < (*edges)[0].size(); j++) {
      all_nodes.push_back((*edges)[i][j]);
    }
  }
  std::set<int> unodes (all_nodes.begin(), all_nodes.end());
  return unodes.size();  
}

void 
addFeasibleEdges (std::vector<std::vector<int> > *edges, 
                  std::vector<int> *deltas, 
                  int epsilon, 
                  int delta,
                  igraph_t *g, 
                  igraph_vector_t *weights, 
                  std::vector<bool> *feasible,
                  std::vector<bool> *already_connected) 
{
  std::vector<bool> newly_feasible (deltas->size(), false);
  int num_newfeasible;
  
  //time(&starttime);
  connected(edges, deltas, g, weights, delta, epsilon, already_connected);
  //time(&endtime);
  //std::cout << "time taken in connected: " << endtime-starttime << std::endl;

  for (int i = 0; i < edges->size(); i++) {
    if ((((*deltas)[i] - epsilon) <= delta) && !((*already_connected)[i]) && !((*feasible)[i])) {
      newly_feasible[i] = true;
      (*feasible)[i] = true;
    }
  }

  num_newfeasible = std::count(newly_feasible.begin(), newly_feasible.end(), true);
  if (num_newfeasible > 0) {
    for (int i = 0; i < edges->size(); i++) {
      if (newly_feasible[i]) {
        igraph_add_edge(g, (*edges)[i][0], (*edges)[i][1]);
        igraph_vector_push_back(weights, (*deltas)[i]);
      }
    }
  }
}

std::set<Combination>
allPairs (int n)
{
  std::set<Combination> combos;
  for (int i = 0; i < (n-1); i++) {
    for (int j = i+1; j < n; j++) {
      combos.insert(Combination(i,j));
    }
  }
  return combos;
}

std::vector<std::vector<int> > 
medianVectors(std::vector<std::vector<int> > triplet) 
{
  if (triplet.size() != 3) { throw(Rcpp::exception("triplet does not have 3 rows.")); }
  int ncols = triplet[0].size();
  
  int index, col, count, cycle, repeat; // for expand.grid into mvecs
  
  std::map<int,int> counts;
  std::map<int,int>::iterator c_it;
  
  std::vector<int> allunique;
  std::map<int,int> majorities;
  
  for (int j = 0; j < ncols; j++) {
    counts.clear();
    for (int i = 0; i < 3; i++) {
      c_it = counts.find(triplet[i][j]);
      if (c_it == counts.end()) {
        counts[triplet[i][j]] = 1; 
      } else {
        counts[triplet[i][j]] += 1;
      }
    }
    if (counts.size() == 3) { 
      allunique.push_back(j); 
    } else {
      for (c_it = counts.begin(); c_it!=counts.end(); ++c_it) {
        if (c_it->second > 1) {
          majorities[j] = c_it->first;
        }
      }
    }
  }
  
  std::vector<std::vector<int> > mvecs (pow(3, allunique.size()), std::vector<int> (ncols));
  cycle = 0;
  for (std::vector<int>::iterator it=allunique.begin(); it != allunique.end(); ++it) {
    repeat = pow(3, cycle);
    count = 0;
    index = 0;
    col = *it;
    for (int i = 0; i < mvecs.size(); i++) {
      if (count == repeat) {
        index += 1;
        count = 0;
      }
      if (index > 2) {
        index = 0;
      }
      mvecs[i][col] = triplet[index][col];
      count += 1;
    }
    cycle += 1;
  }
  for (std::map<int,int>::iterator it=majorities.begin(); it!=majorities.end(); ++it) {
    for (int i = 0; i < mvecs.size(); i++) {
      mvecs[i][it->first] = it->second;
    }
  }
  
  return mvecs;
}


void 
eraseDuplicates (std::vector<std::vector<int> > *new_vectors, 
                 std::vector<std::vector<int> > *existing) 
{
  bool all_equal;
  std::set<int, std::greater<int> > toerase;
  
  for (int k = 0; k < new_vectors->size(); k++) {
    for (int i = 0; i < existing->size(); i++) {
      all_equal = true;
      for (int j = 0; j < (*existing)[i].size(); j++) {
        all_equal = all_equal && (*existing)[i][j] == (*new_vectors)[k][j];
      }
      if (all_equal) {
        toerase.insert(k);
      }
    }
  }
  for (std::set<int>::iterator it=toerase.begin(); it!=toerase.end(); ++it) {
    new_vectors->erase(new_vectors->begin()+(*it));
  }                        
}

std::vector<int> 
calcConnectionCosts(std::vector<std::vector<int> > triplet, 
                    std::vector<std::vector<int> > new_vectors,
                    std::vector<int> *character_weights)
{
  int cost;
  std::vector<int> result;
  
  // iterate through rows of new_vectors
  for (std::vector<std::vector<int> >::iterator xit=new_vectors.begin(); xit!=new_vectors.end(); ++xit) {
    cost = 0;
    // for each row of new_vectors, iterate through all rows of triplet
    for (std::vector<std::vector<int> >::iterator tit=triplet.begin(); tit!=triplet.end(); ++tit) {
      // iterate through columns of both triplets and new_vectors
      for (int i = 0; i < tit->size(); i++) {
        if ((*tit)[i] != (*xit)[i]) {
          cost = cost + (*character_weights)[i];
        }
      }
    }
    result.push_back(cost);
  }
  
  return result;
}

std::vector<std::vector<int> > 
convert2Matrix (IntegerMatrix *x) 
{
  std::vector<std::vector<int> > y (x->nrow(), std::vector<int> (x->ncol()));
  
  for (int i = 0; i < x->nrow(); i++) {
    for (int j = 0; j < x->ncol(); j++) {
      y[i][j] = (*x)(i,j);
    }
  }
  
  return y;
}

std::vector<int> 
convert2Vector (IntegerVector *x)
{
  std::vector<int> y (x->size());
  
  for (int i = 0; i < x->size(); i++) {
    y[i] = (*x)[i];
  }
  
  return y;
}

int 
newSequenceTypes (std::vector<std::vector<int> > *data,
                  std::vector<int> *character_weights,
                  std::vector<std::vector<int> > *edges,
                  std::vector<bool> feasible,
                  int epsilon,
                  std::vector<std::vector<int> > *mvecs)
{
  int row;
  int lambda;
  std::vector<int> v0, v1;
  std::set<Combination> combos;
  std::set<Combination>::iterator cit, pit;
  std::set<int> v_set;
  std::set<std::set<int> > seen_combos;
  std::pair<std::set<std::set<int> >::iterator, bool> ret;
  std::set<int>::iterator v_it;
  std::set<int, std::greater<int> > toerase;

  int ncol = (*data)[0].size();
  int num_starting_mvecs = mvecs->size();

  std::vector<std::vector<int> > new_mvecs;
  std::vector<std::vector<int> > triplet (3, std::vector<int> (ncol));
  
  std::vector<int> costs, new_costs;
  
  for (int i = 0; i < edges->size(); i++) {
    if (feasible[i]) {
      v0.push_back((*edges)[i][0]);
      v1.push_back((*edges)[i][1]);
    }
  }
  
  combos = allPairs(v0.size());
  
  for (cit=combos.begin(); cit!=combos.end(); ++cit) {
    v_set.clear();
    v_set.insert(v0[cit->first]);
    v_set.insert(v0[cit->second]);
    v_set.insert(v1[cit->first]);
    v_set.insert(v1[cit->second]);
    ret = seen_combos.insert(v_set);
    if (ret.second == true && v_set.size() == 3) {
      row = 0;
      for (v_it=v_set.begin(); v_it!=v_set.end(); ++v_it) {
        for (int j = 0; j < ncol; j++) {
          triplet[row][j] = (*data)[*v_it][j];
        }
        row++;
      }
      new_mvecs = medianVectors(triplet);
      eraseDuplicates(&new_mvecs, data);
      eraseDuplicates(&new_mvecs, mvecs);
      new_costs = calcConnectionCosts(triplet, new_mvecs, character_weights);
      costs.insert(costs.end(), new_costs.begin(), new_costs.end());
      mvecs->insert(mvecs->end(), new_mvecs.begin(), new_mvecs.end());
    }
  }
  
  if (costs.size() > 0) {
    lambda = *std::min_element(costs.begin(), costs.end());
    for (int i = 0; i < costs.size(); i++) {
      if (costs[i] > (lambda + epsilon)) {
        toerase.insert(i);
      }
    }
  }

  for (std::set<int>::iterator it=toerase.begin(); it!=toerase.end(); ++it) {
    mvecs->erase(mvecs->begin()+(*it));
  }  
  
  return mvecs->size() - num_starting_mvecs;
}

std::vector<bool> 
deltaStepComponents (std::vector<std::vector<int> > *edges,
                     std::vector<int> *deltas,
                     std::set<int> *unique_deltas,
                     int epsilon) 
{
  igraph_bool_t isconnected;
  int delta, deltaplus;

  std::vector<bool> feasible (edges->size(), false);
  std::vector<bool> already_connected (edges->size(), false);
  
  igraph_t g;
  igraph_empty(&g, nNodes(edges), 0);
  
  igraph_vector_t weights;
  igraph_vector_init(&weights, 0);
  
  for (std::set<int>::iterator it=unique_deltas->begin(); it!=unique_deltas->end(); ++it) {
    delta = *it;
    addFeasibleEdges(edges, deltas, epsilon, delta, &g, &weights, &feasible, &already_connected);
    
    igraph_is_connected(&g, &isconnected, IGRAPH_WEAK);
    if (isconnected) { break; }
  }
  
  if (epsilon > 0) {
    for (std::set<int>::iterator it=unique_deltas->begin(); it!=unique_deltas->end(); ++it) {
      deltaplus = *it;
      if ((deltaplus > delta) && (deltaplus <= (delta+epsilon))) {
        addFeasibleEdges(edges, deltas, epsilon, deltaplus, &g, &weights, &feasible, &already_connected);
      }
    }
  }
    
  igraph_vector_destroy(&weights);
  igraph_destroy(&g);
  
  return feasible;
}

std::set<int, std::greater<int> > 
identifyObsoletes (std::vector<std::vector<int> > *edges,
                   std::vector<bool> *feasible,
                   int num_sampled,
                   int num_total)
{
  std::set<int, std::greater<int> > result;
  std::map<int, int> count;
  
  for (int i = num_sampled; i < num_total; i++) {
    count[i] = 0;
  }
  
  for (int i = 0; i < feasible->size(); i++) {
    if ((*feasible)[i]) {
      if ((*edges)[i][0] >= num_sampled) {
        count[(*edges)[i][0]]++;
      }
      if ((*edges)[i][1] >= num_sampled) {
        count[(*edges)[i][1]]++;
      }
    }
  }

  for (int i = num_sampled; i < num_total; i++) {
    if (count[i] <= 2) {
      result.insert(i);
    }
  }

  return result;
}

void 
mergeAndConvertData(std::vector<std::vector<int> > *alldata, 
                    std::vector<std::vector<int> > *data, 
                    std::vector<int> *character_weights,
                    std::vector<std::vector<int> > *mvecs, 
                    std::vector<int> *deltas, 
                    std::set<int> *unique_deltas,
                    std::vector<std::vector<int> > *edges)
{
  std::vector<std::vector<int> > distances;

  // -- first, merge sampled data with median vectors
  alldata->clear();
  for (std::vector<std::vector<int> >::iterator it=data->begin(); it!=data->end(); ++it) {
    alldata->push_back(*it);
  }
  for (std::vector<std::vector<int> >::iterator it=mvecs->begin(); it!=mvecs->end(); ++it) {
    alldata->push_back(*it);
  }
  // printDataVectors(*alldata);
  // std::cout << std::endl;
  
  // Step 1: Determine the distance matrix for the current
  // sequence types, pool identical sequence types, and order
  // the different distance values as d1 < d2 < ... < dk.

  // (a) determine distance matrix
  distances = hammingDist(alldata, character_weights);
  
  // (b) identical sequence types already pooled (in R)
  
  // -- convert distances matrix to edge list and deltas vector
  deltas->clear();
  edges->clear();
  for (int i = 0; i < (distances.size() - 1); i++) {
    for (int j = i+1; j < distances.size(); j++) {
      deltas->push_back(distances[i][j]);
      std::vector<int> edge;
      edge.push_back(i);
      edge.push_back(j);
      edges->push_back(edge);
    } 
  }
  
  // (c) order unique distance values
  unique_deltas->clear();
  for (std::vector<int>::iterator it=deltas->begin(); it!=deltas->end(); ++it) {
    unique_deltas->insert(*it);
  }
}

void 
mjnC (std::vector<std::vector<int> > *data, 
      std::vector<int> *character_weights,
      std::vector<std::vector<int> > *edgeList, 
      std::vector<int> *edgeLengths, 
      int epsilon)
{
  std::vector<std::vector<int> > alldata;
  std::vector<std::vector<int> > edges; 
  std::vector<int> deltas; 
  std::set<int> unique_deltas;
  std::vector<bool> feasible;
  std::vector<std::vector<int> > mvecs;
  std::set<int, std::greater<int> > obsolete_mvecs;
  
  int num_sampled = data->size();
  int num_new_mvecs = 1;
  int num_obsolete;
  
  // ALGORITHM from Bandelt et al. 1999, pgs 39-40
  
  // PHASE I: Successive selection of median vectors
    
  while (num_new_mvecs > 0) {

    mergeAndConvertData(&alldata, data, character_weights, &mvecs, &deltas, 
                        &unique_deltas, &edges);
    
    // Determine the links between sequence types
    // which describe the (epsilon-relaxed) minimum spanning
    // network.
    //time(&starttime);
    feasible = deltaStepComponents(&edges, &deltas, &unique_deltas, epsilon);
    //time(&endtime);
    //std::cout << "time taken in deltaStepComponents: " << endtime-starttime << std::endl;
    
    // Generate median vectors
    //time(&starttime);
    num_new_mvecs = newSequenceTypes (&alldata, character_weights, &edges, 
                                      feasible, epsilon, &mvecs); 
    //time(&endtime);
    //std::cout << "time taken in newSequenceTypes: " << endtime-starttime << std::endl;
    
    // std::cout << "nnew: " << num_new_mvecs << std::endl;
  } // END PHASE I of ALGORITHM
  
  // PHASE II: Construction of the final network.
  
  // Step 5: Calculate the minimum spanning network for
  // the new set of current sequence types. 
  num_obsolete = 1;
  while (num_obsolete > 0) {
    // identify feasible connections for epsilon = 0
    feasible = deltaStepComponents(&edges, &deltas, &unique_deltas, 0);
    obsolete_mvecs = identifyObsoletes(&edges, &feasible, num_sampled, alldata.size());
    for (std::set<int, std::greater<int> >::iterator it=obsolete_mvecs.begin(); it!=obsolete_mvecs.end(); ++it) {
      mvecs.erase(mvecs.begin()+((*it)-num_sampled));
    }
    mergeAndConvertData(&alldata, data, character_weights, &mvecs, &deltas, &unique_deltas, &edges);
    num_obsolete = obsolete_mvecs.size();
  }
  
  for (int i = 0; i < feasible.size(); i++) {
    if (feasible[i]) {
      edgeList->push_back(edges[i]);
      edgeLengths->push_back(deltas[i]);
    }
  }
}

// [[Rcpp::export]]
List 
mjnRcpp (IntegerMatrix dataR, IntegerVector characterWeights, int epsilon)
{
  List resultList;
  
  std::vector<std::vector<int> > edgeList; 
  std::vector<int> edgeLengths;
  
  std::vector<std::vector<int> > data;
  data = convert2Matrix(&dataR);
  
  std::vector<int> character_weights;
  character_weights = convert2Vector(&characterWeights);
  
  mjnC(&data, &character_weights, &edgeList, &edgeLengths, epsilon);

  IntegerMatrix eList (edgeList.size(), 2);
  IntegerVector eLengths (edgeLengths.size());
  
  for (int i = 0; i < edgeList.size(); i++) {
    eLengths[i] = edgeLengths[i];
    for (int j = 0; j < 2; j++) {
      eList(i,j) = edgeList[i][j];
    }
  }

  resultList["edgeList"] = eList;
  resultList["edgeLengths"] = eLengths;
  return resultList;
}
  