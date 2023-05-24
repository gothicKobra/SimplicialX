#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <set>
#include <vector>
#include <map>
#include <algorithm>
#include <queue>
#include <chrono>
#include <list>
#include <stack>
#include "SimplexInterface.h"
using namespace std;
using namespace arma;
const int TrieMx = 30;
const int TrieINF = 10000000;

class SimplexTrie : public SimpInterface
{
	friend class SimplicialComplex;
private:
	struct TrieNode
	{
		int value;
		int depth;
		int number;
		map<int, TrieNode*> children;
		TrieNode* parent;
		double weight;
	};
	TrieNode* root;
	int f[TrieMx];
	vector<vector<list<TrieNode*>>> depth_lists;
	vector<TrieNode*> nums;
	vector<int>translation;
	SimplexTrie() {
		root = new TrieNode;
		root->depth = -1;
		depth_lists.resize(15000, vector<list<TrieNode*>>(TrieMx));
		for (int i = 0; i < TrieMx; i++) {
			f[i] = 0;
		}
	}
public:
	~SimplexTrie() {
		for (auto c : root->children) {
			delete c.second;
		}
		delete root;
	}

	void addToList(TrieNode* v) {
		depth_lists[v->value][v->depth].push_back(v);
	}

	void deleteFromList(TrieNode* v) {
		int pos = 0;
		auto cur = depth_lists[v->value][v->depth].begin();
		while (*cur != v) {
			cur++;
		}
		depth_lists[v->value][v->depth].erase(cur);
	}

	void dfs(TrieNode* v) {
		for (auto c : v->children) {
			dfs(c.second);
		}
		deleteFromList(v);
		f[v->depth]--;
		delete v;
	}

	bool goodPath(TrieNode* v, vector<int>& word, int pos) {
		if (v->value == word[pos]) {
			pos--;
			if (pos == -1) {
				return true;
			}
		}
		if (v->parent == root) {
			return false;
		}
		return goodPath(v->parent, word, pos);
	}

	void numeration(TrieNode* v, int& t) {
		v->number = t++;
		nums[v->number] = v;
		for (auto c : v->children) {
			numeration(c.second, t);
		}
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

	void bfsFinder(TrieNode* v, vector<int>& word, vector<TrieNode*>& ans) {
		queue<pair<TrieNode*, int>>q;
		q.push({ root, 0 });
		int n = word.size();
		while (!q.empty()) {
			TrieNode* v = q.front().first;
			int pos = q.front().second;
			q.pop();
			if (pos < word.size() && word[pos] == v->value) {
				pos++;
			}
			if (pos == word.size()) {
				ans.push_back(v);
			}
			for (auto c : v->children) {
				q.push({ c.second,pos });
			}
		}
	}

	void openStar(vector<int>& task) {
		vector<int> word = task;
		makePretty(word);
		TrieNode* v = simplexFinder(word);
		if (v == nullptr) {
			throw std::runtime_error("This simplex doesn't exist\n");
		}
		vector<TrieNode*> ans;
		bfsFinder(root, word, ans);
		for (auto c : ans) {
			printSimplex(c);
		}
	}

	void closeDFS(TrieNode* v, set<int>& kids) {
		if (v->children.empty()) {
			TrieNode* cur = v;
			while (cur != root) {
				kids.insert(cur->value);
				cur = cur->parent;
			}
			return;
		}
		for (auto c : v->children) {
			closeDFS(c.second, kids);
		}
	}

	void closeStar(vector<int>& task) {
		vector<int> word = task;
		makePretty(word);
		int n = word.size();
		TrieNode* v = simplexFinder(word);
		if (v == nullptr) {
			throw std::runtime_error("This simplex doesn't exist\n");
		}
		int last = word[n - 1];
		set<int> kids;
		for (int i = n - 1; i < TrieMx; i++) {
			if (depth_lists[last][i].empty()) {
				continue;
			}
			auto now = depth_lists[last][i].begin();
			while (now != depth_lists[last][i].end()) {
				if (!goodPath(*now, word, n - 1)) {
					now++;
					continue;
				}
				TrieNode* temp = *now;
				now++;
				TrieNode* v = temp;
				closeDFS(v, kids);
			}
		}
		queue<TrieNode*>q;
		for (auto c : kids) {
			q.push(root->children[c]);
		}
		vector<TrieNode*> ans;
		while (!q.empty()) {
			TrieNode* v = q.front();
			q.pop();
			ans.push_back(v);
			for (auto c : v->children) {
				if (kids.count(c.first)) {
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
		TrieNode* v = simplexFinder(word);
		if (v == nullptr) {
			throw std::runtime_error("This simplex doesn't exist\n");
		}
		vector<TrieNode*> op;
		bfsFinder(root, word, op);
		int n = word.size();
		int last = word[n - 1];
		set<int> kids;
		for (int i = n - 1; i < TrieMx; i++) {
			if (depth_lists[last][i].empty()) {
				continue;
			}
			auto now = depth_lists[last][i].begin();
			while (now != depth_lists[last][i].end()) {
				if (!goodPath(*now, word, n - 1)) {
					now++;
					continue;
				}
				TrieNode* temp = *now;
				now++;
				TrieNode* v = temp;
				closeDFS(v, kids);
			}
		}
		queue<TrieNode*>q;
		for (auto c : kids) {
			q.push(root->children[c]);
		}
		vector<TrieNode*> cl;
		while (!q.empty()) {
			TrieNode* v = q.front();
			q.pop();
			cl.push_back(v);
			for (auto c : v->children) {
				if (kids.count(c.first)) {
					q.push(c.second);
				}
			}
		}
		int i = 0, j = 0;
		n = op.size();
		int m = cl.size();
		vector<TrieNode*> ans;
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

	void erase(vector<int>& task) {
		vector<int> word = task;
		makePretty(word);
		int n = word.size();
		int last = word[n - 1];
		for (int i = n - 1; i < TrieMx; i++) {
			if (depth_lists[last][i].empty()) {
				continue;
			}
			auto now = depth_lists[last][i].begin();
			while (now != depth_lists[last][i].end()) {
				if (!goodPath(*now, word, n - 1)) {
					now++;
					continue;
				}
				TrieNode* temp = *now;
				now++;
				TrieNode* v = temp;
				v->parent->children.erase(v->value);
				dfs(v);
			}
		}
	}

	void changeWeight(vector<int>& task, double weight = 1) {
		vector<int> word = task;
		makePretty(word);
		TrieNode* ans = simplexFinder(word);
		if (ans == nullptr) {
			std::cout << "This simplex doesn't exist\n";
			return;
		}
		ans->weight = weight;
		cout << "OK\n";
	}

	void insertRecursive(vector<int>& word, TrieNode* v, int pos, double weight) {
		int n = word.size();
		for (int i = pos; i < n; i++) {
			if (!v->children.count(word[i])) {
				TrieNode* temp = new TrieNode;
				v->children[word[i]] = temp;
				temp->value = word[i];
				temp->parent = v;
				temp->depth = v->depth + 1;
				f[v->depth + 1]++;
				temp->weight = 1;
				if (temp->depth == n - 1) {
					temp->weight = weight;
				}
				addToList(temp);
			}
			if (i + 1 < n) {
				insertRecursive(word, v->children[word[i]], i + 1, weight);
			}
		}
	}

	void insert(vector<int>& task, double weight = 1) {
		vector<int> word = task;
		makePretty(word);
		insertRecursive(word, root, 0, weight);
	}

	vector<TrieNode*> simplexList() {
		queue<TrieNode*> q;
		vector<TrieNode*> ans;
		if (root->children.size() == 0) {
			return ans;
		}
		for (auto c : root->children) {
			q.push(c.second);
		}
		while (q.size()) {
			TrieNode* v = q.front();
			q.pop();
			ans.push_back(v);
			for (auto c : v->children) {
				q.push(c.second);
			}
		}
		return ans;
	}

	void fVector() {
		cout << "Current f vector: ";
		for (int i = 0; i < TrieMx; i++) {
			if (f[i] == 0)
				break;
			cout << f[i] << ' ';
		}
		cout << '\n';
	}

	int simplexCount() {
		int ans = 0;
		for (int i = 0; i < TrieMx; i++) {
			if (f[i] == 0)
				break;
			ans += f[i];
		}
		return ans;
	}

	TrieNode* simplexFinder(vector<int>& word) {
		TrieNode* v = root;
		int pos = 0;
		while (true) {
			if (v->value == word[pos]) {
				pos++;
				if (word.size() == pos) {
					return v;
				}
			}
			if (v->children.count(word[pos]) == false) {
				return nullptr;
			}
			v = v->children[word[pos]];
		}
	}

	void dimensionFinder(TrieNode* v, int d, vector<int>& all) {
		if (v->depth == d) {
			all.push_back(v->number);
			return;
		}
		for (auto c : v->children) {
			dimensionFinder(c.second, d, all);
		}
	}

	void dimensionFinderNodes(TrieNode* v, int d, vector<TrieNode*>& all) {
		if (v->depth == d) {
			all.push_back(v);
			return;
		}
		for (auto c : v->children) {
			dimensionFinderNodes(c.second, d, all);
		}
	}
	int bitsCount(int x) {
		int ans = 0;
		for (int i = 1; i < (1 << 30); i *= 2) {
			if (x & i) {
				ans++;
			}
		}
		return ans;
	}

	void findingNeigbours(TrieNode* start, int neighbourDimension, vector<TrieNode*>& all) {
		if (start->depth < neighbourDimension) {
			vector<int> word;
			TrieNode* cur = start;
			while (cur != root) {
				word.push_back(cur->value);
				cur = cur->parent;
			}
			sort(begin(word), end(word));
			int last = start->value;
			for (int i = start->depth; i <= neighbourDimension; i++) {
				if (depth_lists[last][i].empty()) {
					continue;
				}
				auto now = depth_lists[last][i].begin();
				while (now != depth_lists[last][i].end()) {
					if (!goodPath(*now, word, word.size() - 1)) {
						now++;
						continue;
					}
					dimensionFinderNodes(*now, neighbourDimension, all);
					now++;
				}
			}
			return;
		}
		if (start->depth == neighbourDimension) {
			all.push_back(start);
			return;
		}
		vector<int> word;
		TrieNode* cur = start;
		while (cur != root) {
			word.push_back(cur->value);
			cur = cur->parent;
		}
		sort(begin(word), end(word));
		vector<int> res;
		for (int i = 1; i < (1 << word.size()); i++) {
			if (bitsCount(i) == neighbourDimension + 1) {
				res.clear();
				for (int y = 0; y < word.size(); y++) {
					if ((1 << y) & i) {
						res.push_back(word[y]);
					}
				}
				all.push_back(simplexFinder(res));
			}
		}
	}

	vector<vector<pair<int, double>>> graphPQBuilder(int p, int q) {
		int t = 0;
		nums.resize(simplexCount() + 1, nullptr);
		numeration(root, t);
		vector<int> all_p_simplex, all_q_simplex;
		dimensionFinder(root, p, all_p_simplex);
		dimensionFinder(root, q, all_q_simplex);
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
		vector<vector<TrieNode*>> q_neighbours(all_p_simplex.size());
		vector<vector<TrieNode*>> p_neighbours(all_q_simplex.size());
		for (int i = 0; i < all_p_simplex.size(); i++) {
			findingNeigbours(nums[all_p_simplex[i]], q, q_neighbours[i]);
		}
		for (int i = 0; i < all_q_simplex.size(); i++) {
			findingNeigbours(nums[all_q_simplex[i]], p, p_neighbours[i]);
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
		vector<pair<double, int>> dist(edges.size(), { TrieINF, 0 });
		dist[s] = { 0, 1 };
		set<std::pair<int, int>> r;
		r.insert(std::make_pair(dist[s].first, s));
		while (!r.empty()) {
			pair<int, int> v = *r.begin();
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
		return dist;
	}

	double dijkstraSpecial(int s, int t, vector<vector<pair<int, double>>>& edges) {
		vector<pair<double, int>> dist(edges.size(), { TrieINF, 0 });
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
		return TrieINF;
	}

	TrieNode* simplexFinderUpgrade(vector<int>& word, int ind) {
		TrieNode* v = root;
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
			if (v->children.count(word[pos]) == false) {
				return nullptr;
			}
			v = v->children[word[pos]];
		}
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
		TrieNode* v = simplexFinder(wordFirst);
		TrieNode* u = simplexFinder(wordSecond);
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
		return dijkstraSpecial(s, t, graph);
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
			TrieNode* v = nums[rl];
			int ans = graph[i].size();
			vector<int> word;
			while (v != root) {
				word.push_back(v->value);
				v = v->parent;
			}
			for (int y = word.size() - 1; y >= 0; y--) {
				cout << word[y] << ' ';
			}
			cout << '\t' << ans << '\n';
		}
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
		#pragma omp parallel for 
		for (int i = 0; i < n; i++) {
			vector<pair<double, int>> dist = dijkstra(i, graph);
			double sum = 0;
			double cur = 0;
			for (auto c : dist) {
				if (c.first == TrieINF) {
					continue;
				}
				sum += c.first;
				cur++;
			}
			double res = (cur - 1) / sum;
			res *= ((cur - 1) / (ans - 1));
			cent[i] = res;
		}
		vector<int> rev(graph.size());
		for (int i = 0; i < translation.size(); i++) {
			if (translation[i] != -1) {
				rev[translation[i]] = i;
			}
		}
		for (int i = 0; i < graph.size(); i++) {
			int rl = rev[i];
			TrieNode* v = nums[rl];
			vector<int> word;
			while (v != root) {
				word.push_back(v->value);
				v = v->parent;
			}
			for (int y = word.size() - 1; y >= 0; y--) {
				cout << word[y] << ' ';
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
				if (path[s].first == TrieINF) {
					continue;
				}
				dist = dijkstra(s, graph);
				for (int t = s + 1; t < graph.size(); t++) {
					if (t == cur) {
						continue;
					}
					if (path[t].first == TrieINF || path[s].first + path[t].first > dist[t].first) {
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
			std::cout << ans / ((graph.size() - 1) * (graph.size() - 2) / 2) << '\n';
		}
		vector<int> rev(graph.size());
		for (int i = 0; i < translation.size(); i++) {
			if (translation[i] != -1) {
				rev[translation[i]] = i;
			}
		}
		for (int i = 0; i < graph.size(); i++) {
			int rl = rev[i];
			TrieNode* v = nums[rl];
			vector<int> word;
			while (v != root) {
				word.push_back(v->value);
				v = v->parent;
			}
			for (int y = word.size() - 1; y >= 0; y--) {
				cout << word[y] << ' ';
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
				edges.insert({ i, e.first });
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
		for (int i = 0; i < TrieMx; i++) {
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

	void printSimplex(TrieNode* v) {
		TrieNode* cur = v;
		vector<int> all;
		while (cur != root) {
			all.push_back(cur->value);
			cur = cur->parent;
		}
		for (int i = all.size() - 1; i >= 0; i--) {
			cout << all[i] << ' ';
		}
		cout << '\t' << v->weight << '\n';
	}

	void allSimplices() {
		vector<TrieNode*> all = simplexList();
		for (auto c : all) {
			printSimplex(c);
		}
	}

	void rec(int n, int cur, int pos, int sum, vector<int>& res, int k) {
		int x = sum;
		x += (1 << cur);
		if (pos + 1 == k) {
			res.push_back(x);
		}
		else {
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
		if (k >= TrieMx || f[k] == 0) {
			throw std::runtime_error("There aro no k-simplices\n");
			return bd;
		}
		if (f[k - p] == 0) {
			throw std::runtime_error("There aro no (k-p)-simplices\n");
			return bd;
		}
		vector<TrieNode*> k_simp, prev_simp;
		dimensionFinderNodes(root, k, k_simp);
		dimensionFinderNodes(root, k - p, prev_simp);
		if (k_simp.size() == 0) {
			//ошибка
			mat A(1, 1, fill::zeros);
			return A;
		}
		if (prev_simp.size() == 0) {
			//ошибка
			mat A(1, 1, fill::zeros);
			return A;
		}
		int k_len = k_simp.size();
		int prev_len = prev_simp.size();
		if (p == 1) {
			return boundaryMatrixNoPer(k_simp, prev_simp, k, orient);
		}
		mat res(prev_len, k_len, fill::zeros);
		vector<int> comb = combinations(k + 1, k + 1 - p);
		for (int i = 0; i < k_len; i++) {
			vector<int> word;
			TrieNode* v = k_simp[i];
			while (v != root) {
				word.push_back(v->value);
				v = v->parent;
			}
			reverse(begin(word), end(word));
			for (int j : comb) {
				vector<int> vec_out, vec_in;
				for (int y = 0; y < k + 1; y++) {
					if ((1 << y) & j) {
						vec_in.push_back(y);
					}
					else {
						vec_out.push_back(y);
					}
				}
				vector<int> vec_full = vec_out;
				for (auto c : vec_in) {
					vec_full.push_back(c);
				}
				double sign = orient * perm_sign(vec_full);
				TrieNode* v = simplexFinderUpgrade(word, j);
				int ind = 0;
				while (prev_simp[ind] != v) {
					ind++;
				}
				res(ind, i) = sign;
			}
		}
		return res;
	}

	mat boundaryMatrixNoPer(vector<TrieNode*>& k_simp, vector<TrieNode*>& prev_simp, int k, double orient) {
		int k_len = k_simp.size();
		int prev_len = prev_simp.size();
		mat res(prev_len, k_len, fill::zeros);
		TrieNode* simp;
		TrieNode* cur;
		for (int i = 0; i < k_len; i++) {
			double sign = orient;
			simp = k_simp[i];
			auto c = simp->children.begin();
			for (int j = 0; j < k + 1; j++) {
				cur = c->second;
				c++;
				int ind = 0;
				while (prev_simp[ind] != cur) {
					ind++;
				}
				res(ind, i) = sign;
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
		vector<TrieNode*> all;
		dimensionFinderNodes(root, k, all);
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
		dim--;
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
			}
			else {
				B = bounds[k];
				for (int i = 0; i < B.n_cols; i++) {
					if (B.col(i).max() == 0) {
						ker++;
					}
				}
			}
			if (k == dim) {
				im = 0;
			}
			else {
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
