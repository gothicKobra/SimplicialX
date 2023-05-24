#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <iomanip>
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <queue>
#include <chrono>
#include <omp.h>
#include "SimplexInterface.h"
using namespace std;
using namespace arma;

const int HasseMx = 30;
const int HasseINF = 10000000;


class HasseDiagram : private SimpInterface
{
	friend class SimplicialComplex;

private:
	struct HasseNode
	{
		int value;
		int depth;
		int number;
		std::map<int, HasseNode*> faces;
		std::map<int, HasseNode*> cofaces;
		double weight;
		~HasseNode() {
			for (auto c : cofaces) {
				delete c.second;
			}
		}
	};
	HasseNode* root;
	int f[HasseMx];
	vector<HasseNode*> nums;
	vector<int>translation;
	HasseDiagram() {
		root = new HasseNode;
		root->depth = -1;
		root->value = -10;
		for (int i = 0; i < HasseMx; i++) {
			f[i] = 0;
		}
	};
public:
	~HasseDiagram() {
		delete root;
	}

	void makePretty(vector<int>& word) {
		std::set<int> f;
		for (auto c : word) {
			f.insert(c);
		}
		word.clear();
		for (auto c : f) {
			word.push_back(c);
		}
	}

	void dfs(HasseNode* v) {
		if (v == nullptr) {
			return;
		}
		for (auto c : v->faces) {
			c.second->cofaces.erase(c.first);
		}
		f[v->depth]--;
	    while(!v->cofaces.empty()) {
			auto c = v->cofaces.begin();
			dfs(c->second);
		}
		delete v;
	}

	void dfsClose(HasseNode* v, set<int>& ans) {
		if (v->cofaces.empty()) {
			for (auto c : v->faces) {
				ans.insert(c.first);
			}
			return;
		}
		for (auto c : v->cofaces) {
			dfsClose(c.second, ans);
		}
	}

	void bfs(HasseNode* v, vector<HasseNode*>& ans) {
		queue<HasseNode*> q;
		ans.push_back(v);
		for (auto c : v->cofaces) {
			q.push(c.second);
		}
		while (!q.empty()) {
			HasseNode* cur = q.front();
			q.pop();
			ans.push_back(cur);
			for (auto c : cur->cofaces) {
				if (c.first < cur->value) {
					continue;
				}
				q.push(c.second);
			}
		}
	}

	void openStar(vector<int>& task) {
		vector<int> word = task;
		makePretty(word);
		HasseNode* v = simplexFinder(word);
		if (v == nullptr) {
			throw std::runtime_error("This simplex doesn't exist\n");
		}
		vector<HasseNode*> ans;
		bfs(v, ans);
		for (auto c : ans) {
			printSimplex(c);
		}
	}

	void closeStar(vector<int>& task) {
		vector<int> word = task;
		makePretty(word);
		HasseNode* v = simplexFinder(word);
		if (v == nullptr) {
			throw std::runtime_error("This simplex doesn't exist\n");
		}
		set<int> kids;
		dfsClose(v, kids);
		queue<HasseNode*>q;
		for (auto c : kids) {
			q.push(root->cofaces[c]);
		}
		vector<HasseNode*> ans;
		while (!q.empty()) {
			HasseNode* v = q.front();
			q.pop();
			ans.push_back(v);
			for (auto c : v->cofaces) {
				if (c.first > v->value && kids.count(c.first)) {
					q.push(c.second);
				}
			}
		}
		for (auto c : ans) {
			printSimplex(c);
		}
	}
	void link(vector<int>& task) {
		vector<int> word = task;
		makePretty(word);
		vector<HasseNode*> op, cl;
		HasseNode* v = simplexFinder(word);
		if (v == nullptr) {
			throw std::runtime_error("This simplex doesn't exist\n");
		}
		bfs(v, op);
		set<int> kids;
		dfsClose(v, kids);
		queue<HasseNode*>q;
		for (auto c : kids) {
			q.push(root->cofaces[c]);
		}
		while (!q.empty()) {
			HasseNode* v = q.front();
			q.pop();
			cl.push_back(v);
			for (auto c : v->cofaces) {
				if (c.first > v->value&& kids.count(c.first)) {
					q.push(c.second);
				}
			}
		}
		int i = 0, j = 0;
		int n = op.size();
		int m = cl.size();
		vector<HasseNode*> ans;
		while (i < n) {
			if (op[i] == cl[j]) {
				i++;
				j++;
				continue;
			}
			ans.push_back(cl[j]);
			j++;
		}
		while (j < m) {
			ans.push_back(cl[j]);
			j++;
		}
		for (auto c : ans) {
			printSimplex(c);
		}
	}
	void numeration(HasseNode* v, int& t) {
		v->number = t++;
		nums[v->number] = v;
		for (auto c : v->cofaces) {
			if (c.first <= v->value) {
				continue;
			}
			numeration(c.second, t);
		}
	}

	void erase(vector<int>& task) {
		vector<int> word = task;
		if (word.empty()) {
			return;
		}
		makePretty(word);
		HasseNode* v = simplexFinder(word);
		dfs(v);
	}

	void changeWeight(vector<int>& task, double weight = 1) {
		vector<int> word = task;
		makePretty(word);
		HasseNode* ans = simplexFinder(word);
		if (ans == nullptr) {
			std::cout << "This simplex doesn't exist\n";
			return;
		}
		ans->weight = weight;
		cout << "OK\n";
	}

	void insert(vector<int>& task, double weight = 1) {
		vector<int> word = task;
		if (word.empty()) {
			return;
		}
		makePretty(word);
		if (simplexFinder(word) != nullptr) {
			cout << "Simplex already exists\n";
			return;
		}
		int n = word.size();
		std::queue<std::pair<HasseNode*, int>> q;
		for (int i = 0; i < n; i++) {
			if (!root->cofaces.count(word[i])) {
				HasseNode* now = new HasseNode;
				now->value = word[i];
				now->faces[word[i]] = root;
				root->cofaces[word[i]] = now;
				now->depth = 0;
				now->weight = 1;
				if (n == 1) {
					now->weight = weight;
				}
				f[0]++;
			}
			q.push({ root->cofaces[word[i]], 1 << i });
		}
		while (!q.empty()) {
			std::pair<HasseNode*, int> v = q.front();
			q.pop();
			for (int i = 0; i < n; i++) {
				if ((1 << i) <= v.second) {
					continue;
				}
				int sum = v.second + (1 << i);
				if (!v.first->cofaces.count(word[i])) {
					HasseNode* now = new HasseNode;
					now->value = word[i];
					now->depth = 1 + v.first->depth;
					f[now->depth]++;
					for (int y = 0; y < n; y++) {
						if ((v.second & (1 << y)) || y == i) {
							HasseNode* f = simplexFinderUpgrade(word, sum - (1 << y));
							now->faces[word[y]] = f;
							f->cofaces[word[y]] = now;
						}
					}
					now->weight = 1;
					if (sum == (1 << n) - 1) {
						now->weight = weight;
					}
				} 
				q.push({ v.first->cofaces[word[i]], v.second + (1 << i) });
			}
		}
	}

	vector<HasseNode*> simplexList() {
		queue<HasseNode*> q;
		vector<HasseNode*> ans;
		if (root->cofaces.size() == 0) {
			return ans;
		}
		for (auto c : root->cofaces) {
			q.push(c.second);
		}
		while (q.size()) {
			HasseNode* v = q.front();
			q.pop();
			ans.push_back(v);
			for (auto c : v->cofaces) {
				if (c.first < v->value) {
					continue;
				}
				q.push(c.second);
			}
		}
		return ans;
	}

	void fVector() {
		cout << "Current f vector: ";
		for (int i = 0; i < HasseMx; i++) {
			if (f[i] == 0)
				break;
			cout << f[i] << ' ';
		}
		cout << '\n';
	}

	int simplexCount() {
		int ans = 0;
		for (int i = 0; i < HasseMx; i++) {
			if (f[i] == 0)
				break;
			ans += f[i];
		}
		return ans;
	}

	HasseNode* simplexFinder(vector<int>& word) {
		HasseNode* v = root;
		int pos = 0;
		while (true) {
			if (v->value == word[pos]) {
				pos++;
				if (word.size() == pos) {
					return v;
				}
			}
			if (v->cofaces.count(word[pos]) == false) {
				return nullptr;
			}
			v = v->cofaces[word[pos]];
		}
	}

	void dimensionFinder(HasseNode* v, int d, vector<HasseNode*>& all) {
		if (v->depth == d) {
			all.push_back(v);
			return;
		}
		for (auto c : v->cofaces) {
			if (c.first <= v->value) {
				continue;
			}
			dimensionFinder(c.second, d, all);
		}
	}

	HasseNode* simplexFinderUpgrade(vector<int>& word, int ind) {
		HasseNode* v = root;
		int pos = 0;
		while (true) {
			if (((1 << pos) & ind) == false) {
				pos++;
				continue;
			}
			if (v->value == word[pos]) {
				pos++;
				if ((1 << (pos) > ind)) {
					return v;
				}
				continue;
			}
			if (v->cofaces.count(word[pos]) == false) {
				return nullptr;
			}
			v = v->cofaces[word[pos]];
		}
	}

	void findingNeigbours(HasseNode* start, int neighbourDimension, vector<HasseNode*>& all,
		                  vector<int>& visited, int step) {
		if (start->depth < neighbourDimension) {
			std::queue<HasseNode*>q;
			q.push(start);
			visited[start->number] = step;
			while (!q.empty()) {
				HasseNode* v = q.front();
				q.pop();
				if (v->depth == neighbourDimension) {
					all.push_back(v);
					continue;
				}
				for(auto c:v->cofaces){
					if (visited[c.second->number] == step) {
						continue;
					}
					visited[c.second->number] = step;
					q.push(c.second);
				}
			}
			return;
		}
		if (start->depth == neighbourDimension) {
			all.push_back(start);
			return;
		}
		std::queue<HasseNode*>q;
		q.push(start);
		visited[start->number] = step;
		while (!q.empty()) {
			HasseNode* v = q.front();
			q.pop();
			if (v->depth == neighbourDimension) {
				all.push_back(v);
				continue;
			}
			for (auto c : v->faces) {
				if (visited[c.second->number] == step) {
					continue;
				}
				visited[c.second->number] = step;
				q.push(c.second);
			}
		}
		return;
	}

	void vertexDegreePQ(int p, int q) {
		if (p < 0 || q < 0) {
			throw std::runtime_error("p and q can't be negative\n");
			return;
		}
		if (p == q) {
			throw std::runtime_error("p and q can't be equal\n");
			return;
		}
		if (f[p] == 0) {
			throw std::runtime_error("There are no p-simplices\n");
			return;
		}
		if (f[q] == 0) {
			throw std::runtime_error("There are no q-simplices\n");
			return;
		}
		vector<vector<pair<int, double>>> graph = graphPQBuilder(p, q);
		vector<int> rev(graph.size());
		for (int i = 0; i < translation.size(); i++) {
			if (translation[i] != -1) {
				rev[translation[i]] = i;
			}
		}
		for (int i = 0; i < graph.size(); i++) {
			int rl = rev[i];
			HasseNode* v = nums[rl];
			int ans = graph[i].size();
			while (v != root) {
				cout << v->faces.begin()->first << ' ';
				v = v->faces.begin()->second;
			}
			cout << '\t' << ans << '\n';
		}
	}

	vector<vector<pair<int, double>>> graphPQBuilder(int p, int q) {
		int t = 0;
		nums.resize(simplexCount() + 1, nullptr);
		numeration(root, t);
		vector<HasseNode*> a, b;
		vector<int> all_p_simplex, all_q_simplex;
		dimensionFinder(root, p, a);
		dimensionFinder(root, q, b);
		for (auto c : a) {
			all_p_simplex.push_back(c->number);
		}
		for (auto c : b) {
			all_q_simplex.push_back(c->number);
		}
		int mx = 0;
		int mxForQ = 0;
		for (auto c : all_p_simplex) {
			mx = max(mx, c);
		}
		for (auto c : all_q_simplex) {
			mxForQ = max(mxForQ, c);
		}
		translation.resize(mx + 1, -1);
		vector<vector<pair<int, double>>> ans(all_p_simplex.size());
		vector<int> Qtrans(mxForQ + 1, -1);
		for (int i = 0; i < all_p_simplex.size(); i++) {
			translation[all_p_simplex[i]] = i;
		}
		for (int i = 0; i < all_q_simplex.size(); i++) {
			Qtrans[all_q_simplex[i]] = i;
		}
		vector<vector<HasseNode*>> q_neighbours(all_p_simplex.size());
		vector<vector<HasseNode*>> p_neighbours(all_q_simplex.size());
		int step = 1;
		std::vector<int> visited(nums.size(), 0);
		for (int i = 0; i < all_p_simplex.size(); i++) {
			findingNeigbours(nums[all_p_simplex[i]], q, q_neighbours[i], visited, step);
			step++;
		}
		for (int i = 0; i < all_q_simplex.size(); i++) {
			findingNeigbours(nums[all_q_simplex[i]], p, p_neighbours[i], visited, step);
			step++;
		}
		vector<int> cnt(all_p_simplex.size(), 0);
		int shag = 0;
		vector<pair<double, int>> res(all_p_simplex.size(), { 0, 0 });
		for (int i = 0; i < all_p_simplex.size(); i++) {
			shag++;
			for (auto c : q_neighbours[i]) {
				int pos = Qtrans[c->number];
				for (auto e : p_neighbours[pos]) {
					if (cnt[translation[e->number]] < shag) {
						cnt[translation[e->number]] = shag;
						int y = e->number;
						if (i != translation[y]) {
							res[translation[y]] = { c->weight, shag };
						}
					}
					if (cnt[translation[e->number]] == shag) {
						int y = e->number;
						if (i != translation[y] && c->weight < res[translation[y]].first) {
							res[translation[y]] = { c->weight, shag };
						}
					}
				}
			}
			for (int y = 0; y < res.size(); y++) {
				if (res[y].second == shag) {
					ans[i].push_back({ y, res[y].first });
				}
			}
		}
		return ans;
	}

	vector<pair<double, int>> dijkstra(int s, vector<vector<pair<int, double>>>& edges) {
		vector<pair<double, int>> dist(edges.size(), { HasseINF, 0 });
		dist[s] = { 0, 1 };
		set<std::pair<double, int>> r;
		r.insert(std::make_pair(dist[s].first, s));
		while (!r.empty()) {
			pair<double, int> v = *r.begin();
			for (auto c : edges[v.second] ) {
				if (dist[c.first].first == c.second + dist[v.second].first) {
					dist[c.first].second += dist[v.second].second;
				}
				if (dist[c.first].first > c.second + dist[v.second].first) {
					r.erase({ dist[c.first].first, c.first });
					dist[c.first].first = c.second + dist[v.second].first;
					dist[c.first].second = dist[v.second].second;
					r.insert({ dist[c.first].first, c.first });
				}
			}
			r.erase(v);
		}
		return dist;
	}

	double dijkstraSpecial(int s, int t, vector<vector<pair<int, double>>>& edges) {
		vector<pair<double, int>> dist(edges.size(), { HasseINF, 0 });
		dist[s] = { 0, 1 };
		set<std::pair<double, int>> r;
		r.insert(std::make_pair(dist[s].first, s));
		while (!r.empty()) {
			pair<double, int> v = *r.begin();
			if (v.second == t) {
				return v.first;
			}
			for (auto c : edges[v.second]) {
				if (dist[c.first].first == c.second + dist[v.second].first) {
					dist[c.first].second += dist[v.second].second;
				}
				if (dist[c.first].first > c.second + dist[v.second].first) {
					r.erase({ dist[c.first].first, c.first });
					dist[c.first].first = c.second + dist[v.second].first;
					dist[c.first].second = dist[v.second].second;
					r.insert({ dist[c.first].first, c.first });
				}
			}
			r.erase(v);
		}
		return HasseINF;
	}

	double distansePQ(vector<int>& taskFirst, vector<int>& taskSecond, int p, int q) {
		if (p < 0 || q < 0) {
			throw std::runtime_error("p and q can't be negative\n");
			return -1;
		}
		if (p == q) {
			throw std::runtime_error("p and q can't be equal\n");
			return -1;
		}
		if (f[p] == 0) {
			throw std::runtime_error("There are no p-simplices\n");
			return -1;
		}
		if (f[q] == 0) {
			throw std::runtime_error("There are no q-simplices\n");
			return -1;
		}
		if (taskFirst.size() != p + 1 || taskSecond.size() != p + 1) {
			throw std::runtime_error("Simplices' dimensions must be p\n");
			return -1;
		}
		vector<int> wordFirst = taskFirst;
		vector<int> wordSecond = taskSecond;
		makePretty(wordFirst);
		makePretty(wordSecond);
		HasseNode* v = simplexFinder(wordFirst);
		HasseNode* u = simplexFinder(wordSecond);
		if (v == nullptr) {
			throw std::runtime_error("Error: first simplex doesn't exist\n");
			return -1;
		}
		if (u == nullptr) {
			throw std::runtime_error("Error: second simplex doesn't exist\n");
			return -1;
		}
		if (v == u) {
			return 0;
		}
		vector<vector<pair<int, double>>> graph = graphPQBuilder(p, q);
		int s = translation[v->number];
		int t = translation[u->number];
		vector<pair<double, int>> dist = dijkstra(s, graph);
		return dist[t].first;
	}

	void closeness(int p, int q) {
		if (p < 0 || q < 0) {
			throw std::runtime_error("p and q can't be negative\n");
			return;
		}
		if (p == q) {
			throw std::runtime_error("p and q can't be equal\n");
			return;
		}
		if (f[p] == 0) {
			throw std::runtime_error("There are no p-simplices\n");
			return;
		}
		if (f[q] == 0) {
			throw std::runtime_error("There are no q-simplices\n");
			return;
		}
		std::chrono::time_point<std::chrono::high_resolution_clock> start, finish;
		vector<vector<pair<int, double>>> graph = graphPQBuilder(p, q);
		vector<double> cent(graph.size());
		double ans = graph.size();
		size_t n = graph.size();
		vector<int> rev(graph.size());
		for (int i = 0; i < translation.size(); i++) {
			if (translation[i] != -1) {
				rev[translation[i]] = i;
			}
		}
        #pragma omp parallel for 
 		for (int i = 0; i < n; i++) {
			vector<pair<double, int>> dist = dijkstra(i, graph);
			double sum = 0;
			double cur = 0;
			for (auto c : dist) {
				if (c.first == HasseINF) {
					continue;
				}
				sum += c.first;
				cur++;
			}
			double res = (cur - 1) / sum;
			res *= ((cur - 1) / (ans - 1));
			cent[i] = res;
		}
		for (int i = 0; i < graph.size(); i++) {
			int rl = rev[i];
			HasseNode* v = nums[rl];
			while (v != root) {
				cout << v->faces.begin()->first << ' ';
				v = v->faces.begin()->second;
			}
			cout << '\t' << cent[i] << '\n';
		}
	}

	void betweenness(int p, int q) {
		if (p < 0 || q < 0) {
			throw std::runtime_error("p and q can't be negative\n");
			return;
		}
		if (p == q) {
			throw std::runtime_error("p and q can't be equal\n");
			return;
		}
		if (f[p] == 0) {
			throw std::runtime_error("There are no p-simplices\n");
			return;
		}
		if (f[q] == 0) {
			throw std::runtime_error("There are no q-simplices\n");
			return;
		}
		vector<vector<pair<int, double>>> graph = graphPQBuilder(p, q);
		vector<double> cent(graph.size());
		#pragma omp parallel for 
		for (int i = 0; i < graph.size(); i++) {
			int cur = i;
			vector<pair<double, int>> path = dijkstra(cur, graph);
			vector<pair<double, int>> dist;
			double ans = 0;
			for (int s = 0; s < graph.size() - 1; s++) {
				if (s == cur) {
					continue;
				}
				if (path[s].first == HasseINF) {
					continue;
				}
				dist = dijkstra(s, graph);
				for (int t = s + 1; t < graph.size(); t++) {
					if (t == cur) {
						continue;
					}
					if (path[t].first == HasseINF || path[s].first + path[t].first > dist[t].first) {
						continue;
					}
					double x = 0;
					double y = dist[t].second;
					if (path[s].first + path[t].first == dist[t].first) {
						double x1 = path[s].second;
						double x2 = path[t].second;
						x += x1 * x2;
					}
					ans += x / y;
				}
			}
			cent[i] = ans / ((graph.size() - 1) * (graph.size() - 2) / 2);
		}
		vector<int> rev(graph.size());
		for (int i = 0; i < translation.size(); i++) {
			if (translation[i] != -1) {
				rev[translation[i]] = i;
			}
		}
		for (int i = 0; i < graph.size(); i++) {
			int rl = rev[i];
			HasseNode* v = nums[rl];
			while (v != root) {
				cout << v->faces.begin()->first << ' ';
				v = v->faces.begin()->second;
			}
			cout << '\t' << cent[i] << '\n';
		}
	}

	double clusterCoeff(int p, int q) {
		if (p < 0 || q < 0) {
			throw std::runtime_error("p and q can't be negative\n");
			return -1;
		}
		if (p == q) {
			throw std::runtime_error("p and q can't be equal\n");
			return -1;
		}
		if (f[p] == 0) {
			throw std::runtime_error("There are no p-simplices\n");
			return -1;
		}
		if (f[q] == 0) {
			throw std::runtime_error("There are no q-simplices\n");
			return -1;
		}
		vector<vector<pair<int, double>>> graph = graphPQBuilder(p, q);
		double ans = 0;
		double znam = 0;
		set<pair<int, int>> edges;
		for (int i = 0; i < graph.size(); i++) {
			long long x = graph[i].size();
			x *= (x - 1);
			x /= 2;
			znam += x;
			for (auto e : graph[i]) {
				if (e.first < i) {
					continue;
				}
				edges.insert({ i, e.first});
			}
		}
		if (znam == 0) {
			return 0;
		}
		for (int i = 0; i < graph.size() - 2; i++) {
			for (int y = i + 1; y < graph.size() - 1; y++) {
				if (edges.count({ i, y }) == false) {
					continue;
				}
				for (int j = y + 1; j < graph.size(); j++) {
					if (edges.count({ i, j }) == false) {
						continue;
					}
					if (edges.count({ y, j }) == false) {
						continue;
					}
					ans += 3;
				}
			}
		}
		return ans / znam;
	}

	int eulerNumber() {
		int ans = 0;
		for (int i = 0; i < HasseMx; i++) {
			if (f[i] == 0)
				break;
			if (i % 2) {
				ans -= f[i];
			}
			else {
				ans += f[i];
			}
		}
		return ans;
	}

	void printSimplex(HasseNode* v) {
		HasseNode* cur = v;
		while (cur != root) {
			cout << cur->faces.begin()->first << ' ';
			cur = cur->faces.begin()->second;
		}
		cout << '\t' << v->weight << '\n';
	}

	void allSimplices() {
		vector<HasseNode*> all = simplexList();
		for (auto c : all) {
			printSimplex(c);
		}

	}

	void rec(int n, int cur, int pos, int sum, vector<int>& res, int k) {
		int x = sum;
		x += (1 << cur);
		if (pos + 1 == k) {
			res.push_back(x);
		} else {
			for (int i = cur + 1; i < n; i++) {
				rec(n, i, pos + 1, x, res, k);
			}
		}
		if (n - cur - 1 + pos < k) {
			return;
		}
		for (int i = cur + 1; i < n; i++) {
			rec(n, i, pos, sum, res, k);
		}
	}

	vector<int> combinations(int n, int k) {
		vector<int> res;
		vector<int> ans;
		rec(n, 0, 0, 0, res, k);
		set<int> f;
		for (auto c : res) {
			if (f.count(c)) {
				continue;
			}
			f.insert(c);
			ans.push_back(c);
		}
		return ans;
	}

	double perm_sign(vector<int>& a) {
		int cnt = 0;
		for (int i = 0; i < a.size() - 1; i++) {
			for (int j = i + 1; j < a.size(); j++) {
				if (a[i] > a[j]) {
					cnt++;
				}
			}
		}
		cnt %= 2;
		if (cnt == 0) {
			return 1;
		}
		return -1;
	}

	mat boundaryMatrix(int k = 1, int p = 1) {
		double orient = 1;
		mat bd(1, 1, fill::zeros);
		if (k < 1) {
			throw std::runtime_error("k must be greater than 0\n");
			return bd;
		}
		if (p < 1) {
			throw std::runtime_error("p must be greater than 0\n");
			return bd;
		}
		if (k - p < 0) {
			throw std::runtime_error("k can't be less than p\n");
			return bd;
		}
		if (k >= HasseMx || f[k] == 0) {
			throw std::runtime_error("There aro no k-simplices\n");
			return bd;
		}
		if (f[k - p] == 0) {
			throw std::runtime_error("There aro no (k-p)-simplices\n");
			return bd;
		}
		vector<HasseNode*> k_simp, prev_simp;
		dimensionFinder(root, k, k_simp);
		dimensionFinder(root, k - p, prev_simp);
		int k_len = k_simp.size();
		int prev_len = prev_simp.size();
		if (p == 1) {
			return boundaryMatrixNoPer(k_simp, prev_simp, k, orient);
		}
		mat res(prev_len, k_len, fill::zeros);
		vector<int> comb = combinations(k + 1, k + 1 - p);
		for (int i = 0; i < k_len; i++) {
			vector<int> word;
			HasseNode* v = k_simp[i];
			while (v != root) {
				word.push_back(v->faces.begin()->first);
				v = v->faces.begin()->second;
			}
			for (int j : comb) {
				vector<int> vec_out, vec_in;
				for (int y = 0; y < k + 1; y++) {
					if ((1 << y) & j) {
						vec_in.push_back(y);
					} else {
						vec_out.push_back(y);
					}
				}
				vector<int> vec_full = vec_out;
				for (auto c : vec_in) {
					vec_full.push_back(c);
				}
				double sign = orient * perm_sign(vec_full);
				HasseNode* v = simplexFinderUpgrade(word, j);
				int ind = 0;
				while (prev_simp[ind] != v) {
					ind++;
				}
				res(ind, i) = sign;
			}
		}
		return res;
	}

	mat boundaryMatrixNoPer(vector<HasseNode*>& k_simp, vector<HasseNode*>& prev_simp, int k, double orient) {
		int k_len = k_simp.size();
		int prev_len = prev_simp.size();
		mat res(prev_len, k_len, fill::zeros);
		HasseNode* simp;
		HasseNode* cur;
		for (int i = 0; i < k_len; i++) {
			double sign = orient;
			simp = k_simp[i];
			auto c = simp->faces.begin();
			for (int j = 0; j < k + 1; j++) {
				cur = c->second;
				c++;
				int ind = 0;
				while (prev_simp[ind] != cur) {
					ind++;
				}
				res(ind,i) = sign;
				sign = -sign;
			}
		}
		return res;
	}

	mat laplacianMatrix(int k, int p = 1, int q = 1) {
		double orient = 1;
		if (k - p == -1) {
			mat B = boundaryMatrix(k + q, q);
			mat res = B * B.t();
			return res;
		}
		mat B1 = boundaryMatrix(k, p);
		mat B2 = boundaryMatrix(k + q, q);
		mat res = B1.t() * B1 + B2 * B2.t();
		return res;
	}

	mat laplacianMatrixWeight(int k, int p = 1, int q = 1) {
		double orient = 1;
		if (k - p == -1) {
			mat B = boundaryMatrix(k + q, q);
			mat Wk_inv = inv(simplexWeights(k));
			mat Wkq = simplexWeights(k + q);
			mat res = Wk_inv * B * Wkq * B.t();
			return res;
		}
		mat B1 = boundaryMatrix(k, p);
		mat B2 = boundaryMatrix(k + q, q);
		mat Wkp_inv = inv(simplexWeights(k - p));
		mat Wk = simplexWeights(k);
		mat Wk_inv = inv(Wk);
		mat Wkq = simplexWeights(k + q);
		mat res = B1.t() * Wkp_inv * B1 * Wk + Wk_inv * B2 * Wkq * B2.t();
		return res;
	}

	mat simplexWeights(int k) {
		vector<HasseNode*> all; 
		dimensionFinder(root, k, all);
		mat res(all.size(), all.size(), fill::zeros);
		for (int i = 0; i < all.size(); i++) {
			res(i, i) = all[i]->weight;
		}
		return res;
	}

	pair<vec, mat> laplacianSpectre(int k, int p = 1, int q = 1, bool weighted = false) {
		vec val;
		mat res;
		mat lap(1, 1, fill::ones);
		if (weighted) {
			lap = laplacianMatrixWeight(k, p, q);
		}
		else {
			lap = laplacianMatrix(k, p, q);
		}
		arma::eig_sym(val, res, lap);
		return { val, res };
	}

	vec bettiNumbers() {
		int dim = 0;
		while (f[dim] != 0) {
			dim++;
		}
		dim --;
		vec res(dim + 1, fill::zeros);
		vector<mat> bounds;
		mat zer(1, 1, fill::zeros);
		bounds.push_back(zer);
		mat B;
		for (int k = 1; k < dim + 1; k++) {
			B = arma::abs(boundaryMatrix(k));
			B = reduce(B);
			bounds.push_back(B);
		}
		for (int k = 0; k < dim + 1; k++) {
			int ker = 0, im = 0;
			if (k == 0) {
				ker = f[0];
			} else {
				B = bounds[k];
				for (int i = 0; i < B.n_cols; i++) {
					if (B.col(i).max() == 0) {
						ker++;
					}
				}
			}
			if (k == dim) {
				im = 0;
			} else {
				B = bounds[k + 1];
				for (int i = 0; i < B.n_rows; i++) {
					if (B.row(i).max() == 0) {
						continue;
					}
					im++;
				}
			}
			res(k) = ker - im;
		}
		return res;
	}

	mat reduce(mat A, int x = 0) {
		if (x >= A.n_rows || x >= A.n_cols) {
			return A;
		}
		for (int i = x; i < A.n_rows; i++) {
			for (int j = x; j < A.n_cols; j++) {
				if (A(i, j) == 1) {
					A.swap_rows(x, i);
					A.swap_cols(x, j);
					rowvec row_x = A.row(x);
					for (int y = x + 1; y < A.n_rows; y++) {
						if (A(y, x) == 1) {
							rowvec row_y = A.row(y);
							for (int k = 0; k < row_y.size(); k++) {
								row_y(k) += row_x(k);
								if (row_y(k) == 2) {
									row_y(k) = 0;
								}
							}
							A.row(y) = row_y;
						}
					}
					colvec col_x = A.col(x);
					for (int y = x + 1; y < A.n_cols; y++) {
						if (A(x, y) == 1) {
							colvec col_y = A.col(y);
							for (int k = 0; k < col_y.size(); k++) {
								col_y(k) += col_x(k);
								if (col_y(k) == 2) {
									col_y(k) = 0;
								}
							}
							A.col(y) = col_y;
						}
					}
					return reduce(A, x + 1);
				}
			}
		}
		return A;
	}
};
