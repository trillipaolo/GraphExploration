#include <fstream>
#include <queue>
#include <set>
#include <chrono>

using namespace std;
using namespace std::chrono;

short int* levelSynchronousSequentialBFS(bool* graph, int nodes);
void plotLevelTable(short int* lev, int nodes);
set<int> nextNodesSet(bool* graph, int node, int nodes);
queue<int> nextNodesQueue(bool* graph, int node, int nodes);
bool* importGraph(const char* path, const bool undirected);
void printGraph(bool* graph);
int charsToInt(char chars[], int length);

int nodes = 0;

int main(int argc, char* argv[]) {
	bool* graph = importGraph("ERgraphSmall.txt", true);

	//printGraph(graph);

	auto start = high_resolution_clock::now();
	short int* lev = levelSynchronousSequentialBFS(graph, nodes);
	auto stop = chrono::high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(stop - start);
	printf("Sequential BFS: %.2d us\n", duration);


	plotLevelTable(lev, nodes);

	return 0;
}

//short int* levelSynchronousSequentialBFS(bool* graph, int nodes) {
//	/*procedure BFS(r:Node)
//	V = C = 0;; N = root frg.Visited, Current, and Next set
//	r.lev = level = 0
//	repeat
//		C = N
//		for Node c in C do.in parallel
//			for Node n in Nbr(c) do.in parallel
//				if n not in V then
//					N = N + n; V = V + n
//					n.lev = level + 1
//		level++
//	until N = 0;*/
//
//	
//	set<int> visited;
//	set<int> current;
//	set<int> next;
//	set<int>::iterator itr;
//	set<int>::iterator itr2;
//
//	int level = 0;
//	short int* lev = new short int[nodes];
//	for (int i = 0; i < nodes; i++) {
//		lev[i] = -1;
//	}
//
//	// root is node 0
//	visited.insert(0);
//	next.insert(0);
//	lev[0] = 0;
//
//	while (!next.empty()) {
//		current.clear();
//		current = next;
//		next.clear();
//
//		for (itr = current.begin(); itr != current.end(); ++itr) {
//
//			set<int> neighbours = nextNodesSet(graph, *itr, nodes);
//			for (itr2 = neighbours.begin(); itr2 != neighbours.end(); ++itr2) {
//				if (!visited.count(*itr2)) {
//					next.insert(*itr2);
//					visited.insert(*itr2);
//					lev[*itr2] = level + 1;
//				}
//			}
//		}
//		level++;
//	}
//
//	return lev;
//}

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


	bool* visited = new bool[nodes];
	for (int i = 0; i < nodes; i++) {
		visited[i] = false;
	}
	queue<int> current;
	queue<int> next;
	queue<int> neighbours;

	int level = 0;
	short int* lev = new short int[nodes];
	for (int i = 0; i < nodes; i++) {
		lev[i] = -1;
	}

	// root is node 0
	visited[0] = true;
	next.push(0);
	lev[0] = 0;

	int curr;
	int neigh;
	while (!next.empty()) {
		// current = queue<int>();		clean a queue
		swap(current,next);

		while (!current.empty()) {
			curr = current.front();

			neighbours = nextNodesQueue(graph, curr, nodes);
			while (!neighbours.empty()) {
				neigh = neighbours.front();

				if (!visited[neigh]) {
					next.push(neigh);
					visited[neigh] = true;
					lev[neigh] = level + 1;
				}

				neighbours.pop();
			}

			current.pop();
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

queue<int> nextNodesQueue(bool* graph, int node, int nodes) {
	long long int ln = (long long int)node;
	long long int lns = (long long int)nodes;
	queue<int> next;
	for (long long int i = 0; i < nodes; i++) {
		if (graph[ln * lns + i]) {
			next.push(i);
		}
	}

	return next;
}

set<int> nextNodesSet(bool* graph, int node, int nodes) {
	long long int ln = (long long int)node;
	long long int lns = (long long int)nodes;
	set<int> next;
	for (long long int i = 0; i < nodes; i++) {
		if (graph[ln * lns + i]) {
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
	long long int ln = (long long int)nodes;

	// Arrive at edges
	inFile.getline(curr, 256, '#');
	inFile.getline(curr, 256, '\n');

	// Initialize graph without edges
	bool* graph = new bool[ln * ln];
	for (int i = 0; i < nodes * nodes; i++) {
		graph[i] = false;
	}

	// Read edges
	long long int sn;
	long long int en;
	while (!inFile.getline(curr, 256, '\t').eof()) {
		sn = charsToInt(curr, 256);

		inFile.getline(curr, 256, '\n');
		en = charsToInt(curr, 256);

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