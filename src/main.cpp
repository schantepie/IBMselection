// #include <Rcpp.h>
// #include <iostream>
// #include <fstream>
// #include <eigen3/Eigen/Dense>
#include "parent.h"
#include "function_sup.h"
#include <random>
#include <stdint.h>
#include <RcppEigen.h>
#include <math.h>


// using namespace Eigen;
using namespace std;
using namespace Rcpp;

// [[Rcpp::export]]
Rcpp::List ibm(const int nbInd, const int length, const double prob_mut,const int loci,
               const Eigen::MatrixXf & mutationMatrix, const Eigen::MatrixXf & selectionMatrix){


const int traits=2;
const int alleles=2;
const int lociXalleles=loci*alleles ;
const int numberVarCovar=((traits*traits)-traits)/2+traits;

/////////////////////////////////////////////////////////////////////////////////////
// Initialisation of Rmat used to generate multivariate normal for mutations ///////////
/////////////////////////////////////////////////////////////////////////////////////

Eigen::MatrixXf mutCholDecomp(mutationMatrix.llt().matrixL());

/////////////////////////////////////////////////////////////////////////
// Inverse of the selection matrix used to estimate finess  ///////////
//////////////////////////////////////////////////////////////////////

Eigen::MatrixXf selectionMatrixInverse=selectionMatrix.inverse();


//////////////////////////////////////////////////////////////////////////////////////
//////////////// Create the base population of parents and offsprings ///////////////
/////////////////////////////////////////////////////////// /////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

std::vector<parent> p;
for(int i = 0; i < nbInd; ++i)  p.push_back(parent(lociXalleles, traits));

std::vector<parent> o;
for(int i = 0; i < nbInd; ++i) o.push_back(parent(lociXalleles, traits));

////////////////////////////////////////////////////////////////
////////////////// Initialise G-matrix collection //////////////
////////////////////////////////////////////////////////////////

Eigen::MatrixXf Gmat(length,traits*traits);


std::vector<int> femaleRepro(nbInd);
for (int k= 0; k < nbInd; ++k) femaleRepro.push_back(k);

std::vector<double> fitness(nbInd);
std::fill(fitness.begin(), fitness.end(),1);

std::random_device rd;
std::mt19937 generator(rd());

///////////START OF CREATING POPULATION GENERATION

for(int i = 0; i < length; ++i) {

  std::discrete_distribution<int> distribution_discrete(fitness.begin(), fitness.end());

  Eigen::MatrixXf offspringGenetValue4sum(nbInd,traits);

  std::vector<int> maleReproducers(nbInd);
  std::vector<int> femaleReproducers(nbInd);

  ///// Find the male & female reproductors

  for (int f=0;f<nbInd;++f) {
    maleReproducers[f]=distribution_discrete(generator);
    femaleReproducers[f]=RandSamplWithoutReplacSimple(fitness,maleReproducers[f]);
  }

  ///reproduction recombination

  for (int re=0;re<nbInd;++re){
    o[re].reproduction_recombinaison(traits,loci,p[maleReproducers[re]],p[femaleReproducers[re]]);
  }

  ///// mutation , get genetic value, swap indivi

  for (int mu=0;mu<nbInd;++mu){
    o[mu].mutation(prob_mut,traits,lociXalleles,mutCholDecomp);
  }

  ///// get offsprings

  for (int next=0;next<nbInd;++next){
    offspringGenetValue4sum.row(next)=o[next].getGeneticValue();
    p[next].swapIndividual(o[next]);
  }

  ///// estimate fitness of parents for the next repoduction with environmental variance=1
  for (int w=0;w<nbInd;++w){
    double fit=p[w].getGeneticValueEnv(traits).transpose()*selectionMatrixInverse*p[w].getGeneticValueEnv(traits);
    //   double fit=p[w].getGeneticValue().transpose()*selectionMatrixInverse*p[w].getGeneticValue();

    fitness[w]=exp(-fit/2);
  }

  ///// estimate Gmatrix
  MatrixXf centered = offspringGenetValue4sum.rowwise() - offspringGenetValue4sum.colwise().mean();
  MatrixXf cov = (centered.adjoint() * centered) / double(offspringGenetValue4sum.rows() - 1);
  Map<RowVectorXf> VectorG(cov.data(), cov.size());
  Gmat.row(i)=VectorG;
}

return Rcpp::List::create(Rcpp::Named("Gmatrix")=Gmat);

}
