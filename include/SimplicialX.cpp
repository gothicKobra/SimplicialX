#include "SimpStructs.h"

class SimplicialComplex
{
public:
	SimplicialComplex(bool type = true) {
		if (type) {
			complex = new HasseDiagram;
		}
		else {
			complex = new SimplexTrie;
		}
	}
	~SimplicialComplex() {
		delete complex;
	}
	void insert(vector<int> word, double weight = 1) {
		complex->insert(word, weight);
	}
	void erase(vector<int> word) {
		complex->erase(word);
	}
	void changeWeight(vector<int> word, double weight = 1) {
		complex->changeWeight(word, weight);
	}
	int simplexCount() {
		return complex->simplexCount();
	}
	int eulerNumber() {
		return complex->eulerNumber();
	}
	void fVector() {
		complex->fVector();
	}
	void allSimplices() {
		complex->allSimplices();
	}
	void vertexDegreePQ(int p, int q) {
		try {
			complex->vertexDegreePQ(p, q);
		}
		catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			throw;
		}
	}
	void closeness(int p, int q){
		try {
			complex->closeness(p, q);
		}
		catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			throw;
		}
	}
	void betweenness(int p, int q) {
		try {
			complex->betweenness(p, q);
		}
		catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			throw;
		}
	}
	double clusterCoeff(int p, int q) {
		try {
			return complex->clusterCoeff(p, q);
		}
		catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			throw;
		}
	}
	double distancePQ(vector<int> a, vector<int> b, int p, int q) {
		try {
			return complex->distansePQ(a, b, p, q);
		}
		catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			throw;
		}
	}
	mat boundaryMatrix(int k = 1, int p = 1) {
		try {
			return complex->boundaryMatrix(k, p);
		}
		catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			throw;
		}
	}
	mat laplacianMatrix(int k, int p = 1, int q = 1) {
		try {
			return complex->laplacianMatrix(k, p, q);
		}
		catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			throw;
		}
	}
	mat laplacianMatrixWeight(int k, int p = 1, int q = 1) {
		try {
			return complex->laplacianMatrixWeight(k, p, q);
		}
		catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			throw;
		}
	}
	pair <vec, mat> laplacianSpectre(int k, int p = 1, int q = 1, bool weight = false) {
		try {
			return complex->laplacianSpectre(k, p, q, weight);
		}
		catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			throw;
		}
	}
	void openStar(vector<int> word) {
		try {
			complex->openStar(word);
		}
		catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			throw;
		}
	}
	void closeStar(vector<int> word) {
		try {
			complex->closeStar(word);
		}
		catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			throw;
		}
	}
	void link(vector<int> word) {
		try {
			complex->link(word);
		}
		catch (const std::exception & e) {
			std::cerr << e.what() << std::endl;
			throw;
		}
	}
	vec bettiNumbers() {
		return complex->bettiNumbers();
	}
	SimpInterface* complex;
};