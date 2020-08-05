/*
AHN
05.08.2020
by Peng He
*/

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include "igraph.h"

igraph_real_t UNIFORM(float upper) {
  std::random_device rand_dev;
  std::mt19937 generator(rand_dev());
  std::uniform_real_distribution<float>  uni_distr(0, upper);
  return uni_distr(generator);
}

float DISTANCE(float x1, float y1, float x2, float y2) {
  return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

bool FILTER(float d, float lambda, float mu) {
  float a = 1 / (1 + exp(-lambda*(d - mu)));
  float b = UNIFORM(1);
  return a < b;
}

int SAMPLE(igraph_vector_t* V) {
  std::vector<igraph_real_t> Prob(igraph_vector_size(V));
  std::fill(Prob.begin(), Prob.end(), 1);
  std::random_device rand_dev;
  std::mt19937 gen(rand_dev());
  std::discrete_distribution<int> dis_distr(Prob.begin(), Prob.end());
  return VECTOR(*V)[dis_distr(gen)];
}

igraph_t ahn_gen(float A, float L, int N, float lambda, float mu, float eta) {
  igraph_t ahn;
  igraph_i_set_attribute_table(&igraph_cattribute_table);
  igraph_matrix_t mat0, mat;
  igraph_matrix_init(&mat0, N, N);
  igraph_matrix_init(&mat, N, N);
  igraph_vector_t X, Y;
  igraph_vector_init(&X, N);
  igraph_vector_init(&Y, N);
  for(int i = 0; i < N; ++i) {
    VECTOR(X)[i] = UNIFORM(L);
    VECTOR(Y)[i] = UNIFORM(A/L);
  }
  igraph_real_t d;
  for(int i = 0; i < N; ++i) {
    for(int j = i + 1; j < N; ++j) {
      d = DISTANCE(VECTOR(X)[i], VECTOR(Y)[i], VECTOR(X)[j], VECTOR(Y)[j]);
      if(d != 0) {
        igraph_matrix_set(&mat0, i, j, 1 / d);
        igraph_matrix_set(&mat0, j, i, 1 / d);
        if(FILTER(d, lambda, mu)) {
          igraph_matrix_set(&mat, i, j, 1 / d);
          igraph_matrix_set(&mat, j, i, 1 / d);
        }
      }
    }
  }
  igraph_weighted_adjacency(&ahn, &mat, IGRAPH_ADJ_UPPER, NULL, 1);
  igraph_vector_t membership;
  igraph_integer_t no;
  igraph_vector_init(&membership, 0);
  igraph_clusters(&ahn, &membership, NULL, &no, IGRAPH_WEAK);
  while(no > 1) {
    igraph_vector_t rows;
    igraph_vector_t cols;
    igraph_vector_init(&rows, 0);
    igraph_vector_init(&cols, 0);
    igraph_vector_t Comps;
    igraph_vector_init_seq(&Comps, 0, no - 1);
    igraph_integer_t rnd_comp = SAMPLE(&Comps);
    igraph_vector_destroy(&Comps);
    for(int p = 0; p < N; ++p) {
      VECTOR(membership)[p] == rnd_comp ? igraph_vector_push_back(&rows, p) : igraph_vector_push_back(&cols, p);
    }
    igraph_matrix_t temp;
    igraph_matrix_init(&temp, igraph_vector_size(&rows), igraph_vector_size(&cols));
    igraph_matrix_select_rows_cols(&mat0, &temp, &rows, &cols);
    long int temp_i, temp_j;
    igraph_matrix_which_max(&temp, &temp_i, &temp_j);
    igraph_matrix_set(&mat, VECTOR(rows)[temp_i], VECTOR(cols)[temp_j], pow(igraph_matrix_e(&temp, temp_i, temp_j), eta));
    igraph_matrix_set(&mat, VECTOR(cols)[temp_j], VECTOR(rows)[temp_i], pow(igraph_matrix_e(&temp, temp_i, temp_j), eta));
    igraph_weighted_adjacency(&ahn, &mat, IGRAPH_ADJ_UNDIRECTED, NULL, 1);
    igraph_clusters(&ahn, &membership, NULL, &no, IGRAPH_WEAK);    
    igraph_matrix_destroy(&temp);
    igraph_vector_destroy(&rows);
    igraph_vector_destroy(&cols);
  }
  igraph_vector_destroy(&membership);
  igraph_matrix_destroy(&mat0);
  igraph_matrix_destroy(&mat);
  igraph_cattribute_VAN_setv(&ahn, "X", &X);
  igraph_cattribute_VAN_setv(&ahn, "Y", &Y);
  igraph_cattribute_GAS_set(&ahn, "AHN", "AnimalHabitatNetwork");
  igraph_vector_destroy(&X);
  igraph_vector_destroy(&Y);
  return ahn;
}

int main(int argc, char const *argv[]) {
  float A = 25;
  float L = 10;
  int N_patch = 50;
  float Lambda = 30;
  float Mu = .001;
  float Eta = 1;
  // generating an animal habitat network
  igraph_t ahn = ahn_gen(A, L, N_patch, Lambda, Mu, Eta);
  // save the network in graphml format
  FILE* file = fopen("ahn.graphml", "w");
  igraph_write_graph_graphml(&ahn, file, 1);
  fclose(file);
  igraph_destroy(&ahn);
  return 0;
}