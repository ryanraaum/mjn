#include <Rcpp.h>
#include <algorithm> 
#include <iostream> 
#include <igraph/igraph.h>
using namespace Rcpp;

typedef std::pair<int, int> Combination;

IntegerMatrix manhattanDist(IntegerMatrix M) 
{
    int sum;
    int nrows = M.nrow();
    int ncols = M.ncol();
    IntegerMatrix D(nrows);
    for (int i = 0; i < (nrows-1); i++) {
      for (int j = i + 1; j < nrows; j++) {
        sum = 0;
        for (int k = 0; k < ncols; k++) {
          sum = sum + abs(M(i,k) - M(j,k));
        }
        D(i,j) = sum;
        D(j,i) = sum;
      }
    }
    return D;
}

std::vector<bool> connected (IntegerMatrix edges, 
                             const igraph_t *graph_p, 
                             const igraph_vector_t *weights_p, 
                             int delta, int epsilon) 
{
  std::vector<bool> result (edges.nrow(), false);
  igraph_matrix_t pathlengths;
  igraph_matrix_init(&pathlengths, 1, 1);
  igraph_vs_t from;
  igraph_vs_t to;
  igraph_bool_t areconnected;
  int nedges;
  
  for (int i = 0; i < edges.nrow(); i++) {
    igraph_st_edge_connectivity(graph_p, &nedges, edges(i,0), edges(i,1));
    igraph_are_connected(graph_p, edges(i,0), edges(i,1), &areconnected);
    
    if (nedges == 0) {
      result[i] = false;
    } else if ((nedges > 0) && areconnected) {
      result[i] = true;
    } else if (nedges > 0) {
      igraph_vs_1(&from, edges(i,0));
      igraph_vs_1(&to, edges(i,1));
      igraph_shortest_paths_dijkstra(graph_p, &pathlengths, from, to, weights_p, IGRAPH_ALL);
      if (MATRIX(pathlengths, 0, 0) < (delta - epsilon)) {
        result[i] = true;
      }
      igraph_vs_destroy(&from);
      igraph_vs_destroy(&to);
    }
  }
  
  igraph_matrix_destroy(&pathlengths);
  return result;
}

int nnodes (IntegerMatrix edges) 
{
  std::vector<int> all_nodes;
  for (int i=0; i < edges.nrow(); i++) {
    for (int j=0; j < edges.ncol(); j++) {
      all_nodes.push_back(edges(i,j));
    }
  }
  std::set<int> unodes (all_nodes.begin(), all_nodes.end());
  return unodes.size();  
}

void addFeasibleEdges (IntegerMatrix edges, 
                         IntegerVector deltas, 
                         int epsilon, 
                         double delta,
                         igraph_t *g_ptr, 
                         igraph_vector_t *weights_ptr, 
                         std::vector<bool> *feasible_ptr) 
{
  std::vector<bool> already_connected;
  std::vector<bool> newly_feasible (deltas.size(), false);
  int num_newfeasible;
    
  already_connected = connected(edges, g_ptr, weights_ptr, delta, epsilon);
  for (int i = 0; i < edges.nrow(); i++) {
    if ((deltas[i] <= delta) && !(already_connected[i]) && !((*feasible_ptr)[i])) {
      newly_feasible[i] = true;
      (*feasible_ptr)[i] = true;
    }
  }
    
  num_newfeasible = std::count(newly_feasible.begin(), newly_feasible.end(), true);
  if (num_newfeasible > 0) {
    for (int i = 0; i < edges.nrow(); i++) {
      if (newly_feasible[i]) {
        igraph_add_edge(g_ptr, edges(i,0), edges(i,1));
        igraph_vector_push_back(weights_ptr, deltas[i]);
      }
    }
  }
}

std::vector<int> identifyObsoletes (igraph_t *g_ptr, std::vector<int> indices) 
{
  std::vector<int> result;
  igraph_vector_t neighbors;
  igraph_vector_init(&neighbors, 1);
  
  for (std::vector<int>::iterator it=indices.begin(); it!=indices.end(); ++it) {
    igraph_neighbors(g_ptr, &neighbors, *it, IGRAPH_ALL);
    if (igraph_vector_size(&neighbors) <= 2) {
      result.push_back(*it);
    }
  }
  return result;
}

std::set<Combination> allPairs (int n)
{
  std::set<Combination> combos;
  for (int i = 0; i < (n-1); i++) {
    for (int j = i+1; j < n; j++) {
      combos.insert(Combination(i,j));
    }
  }
  return combos;
}

std::set<Combination> convert2Set(IntegerMatrix prev_combos)
{
  std::set<Combination> prior_combos;
  
  if (prev_combos.nrow() > 0 && prev_combos.ncol() != 2) {
    throw(Rcpp::exception("prev_combos does not have 2 columns"));
  }
  for (int i = 0; i < prev_combos.nrow(); i++) {
    prior_combos.insert(Combination(prev_combos(i,0), prev_combos(i,1)));
  }
  
  return prior_combos;
}

std::vector<std::vector<int> > medianVectors(std::vector<std::vector<int> > triplet) 
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

// [[Rcpp::export]]
List newSequenceTypes (IntegerMatrix data, 
                       IntegerMatrix edges, 
                       std::vector<bool> feasible,
                       int epsilon,
                       double lambda,
                       IntegerMatrix prev_combos) 
{
  List returnList;
  int row;
  bool all_equal;
  std::vector<int> v0, v1;
  std::set<Combination> combos, prior_combos;
  std::set<Combination>::iterator it, pit;
  std::set<int> v_set;
  std::set<int>::iterator v_it;
  std::set<int, std::greater<int> > toerase;
  
  std::vector<std::vector<int> > mvecs;
  std::vector<std::vector<int> > triplet (3, std::vector<int> (data.ncol()));
  
  // probably will be refactored out at some point
  prior_combos = convert2Set(prev_combos);

  for (int i = 0; i < edges.nrow(); i++) {
    if (feasible[i]) {
      v0.push_back(edges(i,0));
      v1.push_back(edges(i,1));
    }
  }
  
  combos = allPairs(v0.size());
  if (prior_combos.size() > 0) {
    for (it=combos.begin(); it!=combos.end(); ++it) {
      pit = prior_combos.find(*it);
      if (pit!=prior_combos.end()) {
        combos.erase(*it);
      }
    }
  }
  
  for (it=combos.begin(); it!=combos.end(); ++it) {
    v_set.clear();
    v_set.insert(v0[it->first]);
    v_set.insert(v0[it->second]);
    v_set.insert(v1[it->first]);
    v_set.insert(v1[it->second]);
    if (v_set.size() == 3) {
      row = 0;
      for (v_it=v_set.begin(); v_it!=v_set.end(); ++v_it) {
        for (int j = 0; j < data.ncol(); j++) {
          triplet[row][j] = data(*v_it,j);
        }
        row++;
      }
      mvecs = medianVectors(triplet);
      toerase.clear();
      for (int k = 0; k < mvecs.size(); k++) {
        for (int i = 0; i < data.nrow(); i++) {
          all_equal = true;
          for (int j = 0; j < data.ncol(); j++) {
            all_equal = all_equal && data(i,j) == mvecs[k][j];
          }
          if (all_equal) {
            toerase.insert(k);
          }
        }
      }
      for (std::set<int>::iterator it=toerase.begin(); it!=toerase.end(); ++it) {
        mvecs.erase(mvecs.begin()+(*it));
      }
//      for (std::vector<std::vector<int> >::iterator mit=mvecs.begin(); mit!=mvecs.end(); ++mit) {
//        for (std::vector<int>::iterator rit=mit->begin(); rit!=mit->end(); ++rit) {
//          std::cout << *rit << " ";
//        }
//        std::cout << std::endl;
//      }
      ### at line 124 in mjn.R
    }
  }
  
  returnList["v0"] = v0;
  returnList["v1"] = v1;
  return returnList;
}

// [[Rcpp::export]]
List deltaStepComponents (IntegerMatrix edges, IntegerVector deltas, int epsilon) 
{
  List returnList;
  igraph_bool_t isconnected;
  int delta, deltaplus;
  
  std::vector<bool> feasible (edges.nrow(), false);
  
  igraph_t g;
  igraph_empty(&g, nnodes(edges), 0);
  
  igraph_vector_t weights;
  igraph_vector_init(&weights, 0);
    
  std::set<int> udeltas (deltas.begin(), deltas.end());                         
  
  for (std::set<int>::iterator it=udeltas.begin(); it!=udeltas.end(); ++it) {
    delta = *it;
    addFeasibleEdges(edges, deltas, epsilon, delta, &g, &weights, &feasible);
    
    igraph_is_connected(&g, &isconnected, IGRAPH_WEAK);
    if (isconnected) { break; }
  }
  
  if (epsilon > 0) {
    for (std::set<int>::iterator it=udeltas.begin(); it!=udeltas.end(); ++it) {
      deltaplus = *it;
      if ((deltaplus > delta) && (deltaplus <= (delta+epsilon))) {
        addFeasibleEdges(edges, deltas, epsilon, deltaplus, &g, &weights, &feasible);
      }
    }
  }
  
  igraph_vector_destroy(&weights);
  igraph_destroy(&g);
  
  returnList["delta"] = delta;
  returnList["feasible"] = feasible;
  returnList["epsilon"] = epsilon;
  return returnList;
}
  