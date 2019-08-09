#include "stdafx.h"

#include <fstream>
#include <set>
#include <chrono>

using namespace std;
using namespace std::chrono;

short int* levelSynchronousSequentialBFS(bool* graph, int nodes);
void plotLevelTable(short int* lev, int nodes);
set<int> nextNodes(bool* graph, int node, int nodes);
bool* importGraph(const char* path, const bool undirected);
void printGraph(bool* graph);
int charsToInt(char chars[], int length);

int nodes = 0;

int main(int argc, char* argv[]) {

	bool* graph = importGraph("ERgraphTest.txt", true);

	printGraph(graph);

	auto start = high_resolution_clock::now();
	short int* lev = levelSynchronousSequentialBFS(graph, nodes);
	auto stop = chrono::high_resolution_clock::now();
	auto duration = duration_cast<microseconds>(stop - start);
	printf("Sequential BFS: %d ms\n", duration);


	plotLevelTable(lev, nodes);

	return 0;
}

short int* levelSynchronousSequentialBFS(bool* graph, int nodes) {
	/*procedure BFS(r:Node)
	V = C = 0;; N = root frg.Visited, Current, and Next set
	r.lev = level = 0
	repeat
		C = N
		for Node c in C do.in parallel
			for Node n in Nbr(c) do.in parallel
				if n not in V then
					N = N + n; V = V + n
					n.lev = level + 1
		level++
	until N = 0;*/

	
	set<int> visited;
	set<int> current;
	set<int> next;
	set<int>::iterator itr;
	set<int>::iterator itr2;

	int level = 0;
	short int* lev = new short int[nodes];
	for (int i = 0; i < nodes; i++) {
		lev[i] = -1;
	}

	// root is node 0
	visited.insert(0);
	next.insert(0);
	lev[0] = 0;

	while (!next.empty()) {
		current.clear();
		current = next;
		next.clear();

		for (itr = current.begin(); itr != current.end(); ++itr) {

			set<int> neighbours = nextNodes(graph, *itr, nodes);
			for (itr2 = neighbours.begin(); itr2 != neighbours.end(); ++itr2) {
				if (!visited.count(*itr2)) {
					next.insert(*itr2);
					visited.insert(*itr2);
					lev[*itr2] = level + 1;
				}
			}
		}
		level++;
	}

	return lev;
}

void plotLevelTable(short int* lev, int nodes) {
	int maxLev = -1;
	for (int i = 0; i < nodes; i++) {
		if (lev[i] > maxLev) {
			maxLev = lev[i];
		}
	}

	int* count = new int[maxLev + 2];
	for (int i = 0; i < maxLev + 2; i++) {
		count[i] = 0;
	}
	for (int i = 0; i < nodes; i++) {
		count[lev[i] + 1]++;
	}

	printf("________________________________\n");
	printf("\nLevel\tNodes\tPercentage\n");
	printf("________________________________\n\n");
	for (int i = 1; i < maxLev + 2; i++) {
		printf("%d\t%d\t%.2f %%\n", i - 1, count[i], (float)((count[i]) * 100) / nodes);
	}
	printf("________________________________\n\n");
	printf("Total:\t%d\t%.2f %%\n", nodes - count[0], (float)((nodes - count[0]) * 100) / nodes);
	printf("________________________________\n\n");
}

set<int> nextNodes(bool* graph, int node, int nodes) {
	set<int> next;
	for (int i = 0; i < nodes; i++) {
		if (graph[node * nodes + i]) {
			next.insert(i);
		}
	}

	return next;
}

bool* importGraph(const char* path, const bool undirected) {
	ifstream inFile;
	inFile.open(path);
	if (!inFile) {
		printf("Unable to open file: %s not found.", path);
		exit(1); // terminate with error
	}

	char curr[256];
	// Arrive at value of Nodes:
	inFile.getline(curr, 256, '#');
	inFile.getline(curr, 256, '#');
	inFile.getline(curr, 256, ':');
	inFile.getline(curr, 256, ' ');

	// Read nodes and save it as global variable
	inFile.getline(curr, 256, ' ');
	nodes = charsToInt(curr, 256);

	// Arrive at edges
	inFile.getline(curr, 256, '#');
	inFile.getline(curr, 256, '\n');

	// Initialize graph without edges
	bool* graph = new bool[nodes*nodes];
	for (int i = 0; i < nodes * nodes; i++) {
		graph[i] = false;
	}

	// Read edges
	while (!inFile.getline(curr, 256, '\t').eof()) {

		int sn = charsToInt(curr, 256);
		inFile.getline(curr, 256, '\n');
		int en = charsToInt(curr, 256);

		graph[sn * nodes + en] = true;
		if (undirected) {
			graph[en * nodes + sn] = true;
		}
	}

	inFile.close();

	return graph;
}

void printGraph(bool* graph) {
	for (int i = 0; i < nodes; i++) {
		for (int j = 0; j < nodes; j++) {
			if (graph[i * nodes + j]) {
				printf("1 ");
			}
			else {
				printf("0 ");
			}
		}
		printf("\n");
	}
}

int charsToInt(char chars[], int length) {
	bool cond = true;
	int dim = 0;
	for (int i = 0; i < length && cond; i++) {
		if (chars[i] == '\0') {
			dim = i;
			cond = false;
		}
	}

	if (cond) {
		printf("There is no terminator.");
		exit(1);
	}

	int output = 0;
	for (int i = dim - 1, digit = 1; i >= 0; i--, digit *= 10) {
		if (chars[i] == '-') {
			return -1 * output;
		}
		output += (chars[i] - '0') * digit;
	}

	return output;
}