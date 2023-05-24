#pragma once
#include <vector>
#include <armadillo>
using namespace std;
class SimplicialComplex;
class Trie;
class HasseDiagram;
class SimpInterface {
public:
    virtual void insert(std::vector<int>& word, double weight = 1) = 0;
    virtual void erase(std::vector<int>& word) = 0;
    virtual void changeWeight(std::vector<int>& word, double weight = 1) = 0;
    virtual int simplexCount() = 0;
    virtual void fVector() = 0;
    virtual void vertexDegreePQ(int p, int q) = 0;
    virtual double distansePQ(std::vector<int>& wordFirst, std::vector<int>& wordSecond, int p, int q) = 0;
    virtual void allSimplices() = 0;
    virtual void closeness(int p, int q) = 0;
    virtual int eulerNumber() = 0;
    virtual void betweenness(int p, int q) = 0;
    virtual double clusterCoeff(int p, int q) = 0;
    virtual arma::mat boundaryMatrix(int k = 1, int p = 1) = 0;
    virtual arma::mat laplacianMatrix(int k, int p = 1, int q = 1) = 0;
    virtual arma::mat laplacianMatrixWeight(int k, int p = 1, int q = 1) = 0;
    virtual pair<arma::vec, arma::mat> laplacianSpectre(int k, int p = 1, int q = 1, bool weighted = false) = 0;
    virtual arma::vec bettiNumbers() = 0;
    virtual void openStar(vector<int>& word) = 0;
    virtual void closeStar(vector<int>& word) = 0;
    virtual void link(vector<int>& word) = 0;
private:
    SimpInterface() {};
    friend class HasseDiagram;
    friend class SimplexTrie;
};