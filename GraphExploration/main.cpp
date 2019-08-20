#include <fstream>
#include <queue>
#include <set>
#include <chrono>
#include <thread>
#include <mutex>
#include <omp.h>

using namespace std;
using namespace std::chrono;


// SHARED MEMORY //

int concurentThreadsSupported;
mutex mtx;
enum state { SEQ, QUEUE, QUEUE_TO_READ, READ, READ_TO_QUEUE, END };

// Variables for state change in Hybrid BFS algorithm
int t1 = 64;
int t2 = 262144;
int t3 = 2048;
float alpha = 2.0;
float beta = 2.0;


// FUNCTIONS IDENTIFIER //

short int* sequentialBasedBFS();

short int* readBasedBFS();
void readBFSarray(int num);
void readBFScurrent(int curr);

short int* queueBasedBFS();
void queueBFSarray(int num, bool* visited, vector<int>& next);
void queueBFScurrent(int curr, bool* visited, vector<int>& next);
void queueBFSneighbours(int neigh, bool* visited, vector<int>& next);

void initializeVisitedLev(int num, bool* visited, short int* lev);
void plotLevelTable(short int* lev, int n_nodes);
queue<int> nextNodesQueue(bool* graph, int node, int n_nodes);
vector<int> nextNodesVector(bool* graph, int node, int nodes);
bool* importGraph(const char* path, const bool undirected, int &n_nodes);
void printGraph(bool* graph, int nodes);
int charsToInt(char chars[], int length);


// FUNCTIONS

short int* hybridBFS() {
	state nextState = SEQ;

	int n_curr = 0;
	int n_next = 0;
	bool exp = false;

	bool cond = true;
	while (cond) {
		switch (nextState) {
		case SEQ:

			if (n_next > t1) {
				nextState = QUEUE;
			}
			break;

		case QUEUE:

			if (n_next > t2 || (n_next > alpha * n_curr)) {
				nextState = QUEUE_TO_READ;
				exp = n_next > alpha * n_curr;
			}
			break;

		case QUEUE_TO_READ:

			nextState = READ;
			break;

		case READ:

			if (n_next < t2 || (exp && (n_next < beta * n_curr))) {
				nextState = READ_TO_QUEUE;
			}
			break;

		case READ_TO_QUEUE:

			if (n_next > t1) {
				nextState = QUEUE;
			}
			else {
				nextState = SEQ;
			}

			break;
		default:
			printf("Default branch!\n");
			break;
		}
	}

	return NULL;
}

short int* sequentialBasedBFS(bool* graph, const int n_nodes, const int root) {
	//procedure BFS(r:Node)
	//V = C = 0;; N = root frg.Visited, Current, and Next set
	//r.lev = level = 0
	//repeat
	//	C = N
	//	for Node c in C do
	//		for Node n in Nbr(c) do
	//			if n not in V then
	//				N = N + n; V = V + n
	//				n.lev = level + 1
	//	level++
	//until N = 0;

	// Initialize visited[] and lev[]
	short int level = 0;
	bool* visited = new bool[n_nodes];
	short int* lev = new short int[n_nodes];

	for (int i = 0; i < n_nodes; i++) {
		visited[i] = false;
		lev[i] = -1;
	}

	queue<int> current;
	queue<int> next;
	queue<int> neighbours;

	// root is node 0
	visited[root] = true;
	next.push(root);
	lev[root] = 0;

	int curr;
	int neigh;
	while (!next.empty()) {
		// current = queue<int>();		clean a queue
		swap(current,next);

		while (!current.empty()) {
			curr = current.front();

			neighbours = nextNodesQueue(graph, n_nodes, curr);
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

	delete[] visited;

	return lev;
}

//short int* readBasedBFS(bool* graph, int n_nodes, int root) {
//	//BFS_Read(G: Graph, r : Node) {
//	//	Bitmap V;
//	//	Bool fin[threads];
//	//	V.set(r.id);
//	//	int level = 0; r.lev = level;
//	//	bool finished = false;
//	//	while (!finished) {
//	//		fork;
//	//		fin[tid] = true;
//	//		foreach(c: G.Nodes.partition(tid)) {
//	//			if (c.lev != level)		// if it is not current, go to the next node
//	//				continue;
//	//			foreach(n: c.nbrs) {
//	//				if (!V.isSet(n.id)) { // test and test-and-set
//	//					if (V.atomicSet(n.id)) {
//	//						n.lev = level + 1;
//	//						fin[tid] = false;
//	//					}
//	//				}
//	//			}
//	//		}
//	//		join;
//	//		finished = logicalAnd(fin, threads);
//	//		level++;
//	//	}
//	//}
//	
//	// Initialize visited[] and lev[]
//	short int level = 0;
//	bool* visited = new bool[n_nodes];
//	short int* lev = new short int[n_nodes];
//
//	//#pragma omp parallel for
//	for (int i = 0; i < n_nodes; i++) {
//		visited[i] = false;
//		lev[i] = -1;
//	}
//
//	// root is node 0
//	visited[root] = true;
//	lev[root] = 0;
//
//	bool cond = true;
//	while (cond) {
//		cond = false;
//		for (int i = 0; i < concurentThreadsSupported; i++) {
//			at[i] = thread(readBFSarray, i);
//		}
//		for (int i = 0; i < concurentThreadsSupported; i++) {
//			at[i].join();
//		}
//
//		shared_level++;
//	}
//
//	delete[] at;
//	delete[] visited;
//
//	return shared_lev;
//}
//
//void readBFSarray(int num) {
//	int start = (shared_nodes * num) / concurentThreadsSupported;
//	int end = (shared_nodes * (num + 1)) / concurentThreadsSupported;
//	if (end > shared_nodes) {
//		end = shared_nodes;
//	}
//
//	queue<thread> threads;
//
//	for (int i = start; i < end; i++) {
//		if (shared_lev[i] == shared_level) {
//			shared_cond = true;
//			threads.push(thread(readBFScurrent, i));
//		}
//	}
//		
//	while (!threads.empty()) {
//		threads.front().join();
//		threads.pop();
//	}
//}
//
//void readBFScurrent(int curr) {
//	queue<int> neighbours;
//	neighbours = nextNodesQueue(shared_graph, curr, shared_nodes);
//
//	int neigh;
//	while (!neighbours.empty()) {
//
//		neigh = neighbours.front();
//
//		if (shared_lev[neigh] == -1) {
//			shared_lev[neigh] = shared_level + 1;
//		}
//
//		neighbours.pop();
//	}
//}

short int* queueBasedBFS(bool* graph, const int n_nodes, const int root) {

	vector<int> current;
	vector<int> next;
	queue<thread> threads;

	// Initialize visited[] and lev[]
	short int level = 0;
	bool* visited = new bool[n_nodes];
	short int* lev = new short int[n_nodes];

	for (int i = 0; i < n_nodes; i++) {
		visited[i] = false;
		lev[i] = -1;
	}

	// root is node 0
	visited[root] = true;
	next.push_back(root);
	lev[root] = 0;

	while (!next.empty()) {

		current.clear();
		current.swap(next);

		#pragma omp parallel for
		for (int i = 0; i < current.size(); i++) {
			int curr = current.at(i);

			vector<int> neighbours = nextNodesVector(graph, n_nodes, curr);

			#pragma omp parallel for
			for (int j = 0; j < neighbours.size(); j++) {
				int neigh = neighbours.at(j);
				if (!visited[neigh]) {
					visited[neigh] = true;
					lev[neigh] = level + 1;

					#pragma omp critical (nextPush)
					{
						next.push_back(neigh);
					}
				}
			}
		}

		level++;
	}

	delete[] visited;

	return lev;
}

//void queueBFScurrent(int curr, bool* visited, vector<int> &next) {
//	vector<int> neighbours;
//	queue<thread> threads;
//
//	int neigh;
//
//	neighbours = nextNodesVector(shared_graph, curr, shared_nodes);
//	while (!neighbours.empty()) {
//
//		neigh = neighbours.back();
//
//		if (!visited[neigh]) {
//			threads.push(thread(queueBFSneighbours, neigh, visited, ref(next)));
//		}
//
//		neighbours.pop_back();
//	}
//
//	while (!threads.empty()) {
//		threads.front().join();
//		threads.pop();
//	}
//}

void queueBFSneighbours(int neigh, bool* visited, vector<int> &next, short int* lev, short int level) {
	if (!visited[neigh]) {
		visited[neigh] = true;
		lev[neigh] = level + 1;

		mtx.lock();
		next.push_back(neigh);
		mtx.unlock();
	}
}

//void initializeVisitedLev(int num, bool* visited, short int* lev) {
//	int start = (shared_nodes * num) / concurentThreadsSupported;
//	int end = (shared_nodes * (num + 1)) / concurentThreadsSupported;
//	if (end > shared_nodes) {
//		end = shared_nodes;
//	}
//	for (int i = start; i < end; i++) {
//		visited[i] = false;
//		lev[i] = -1;
//	}
//}

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

	delete[] count;
}

queue<int> nextNodesQueue(bool* graph, int n_nodes, int curr) {
	// TODO can be parallelized

	long long int lc = (long long int)curr;
	long long int ln_ns = (long long int)n_nodes;
	queue<int> next;
	for (long long int i = 0; i < ln_ns; i++) {
		if (graph[lc * ln_ns + i]) {
			next.push(i);
		}
	}

	return next;
}

vector<int> nextNodesVector(bool* graph, int n_nodes, int curr) {

	long long int lc = (long long int)curr;
	long long int ln_ns = (long long int)n_nodes;
	vector<int> neighbours;

	#pragma omp parallel for
	for (long long int i = 0; i < ln_ns; i++) {
		if (graph[lc * ln_ns + i]) {
			#pragma omp critical (neighboursPush)
			{
				neighbours.push_back(i);
			}
		}
	}

	return neighbours;
}

bool* importGraph(const char* path, const bool undirected, int* n_nodes) {
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
	*n_nodes = charsToInt(curr, 256);
	long long int ln = (long long int)*n_nodes;

	// Arrive at edges
	inFile.getline(curr, 256, '#');
	inFile.getline(curr, 256, '\n');

	// Initialize graph without edges
	bool* graph = new bool[ln * ln];
	for (long long int i = 0; i < ln * ln; i++) {
		graph[i] = false;
	}

	// Read edges
	long long int sn;
	long long int en;
	while (!inFile.getline(curr, 256, '\t').eof()) {
		sn = charsToInt(curr, 256);

		inFile.getline(curr, 256, '\n');
		en = charsToInt(curr, 256);

		graph[sn * ln + en] = true;
		if (undirected) {
			graph[en * ln + sn] = true;
		}
	}

	inFile.close();

	return graph;
}

void printGraph(bool* graph, int nodes) {
	printf("\nGraph:\n");
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
	printf("\n");
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


// MAIN

int main(int argc, char* argv[]) {

	// For time computation
	time_point<steady_clock> start;
	time_point<steady_clock> stop;
	milliseconds duration;


	// Import the graph from file
	const int root = 0;
	short int* lev;
	int n_nodes = -1;
	printf("Start Importing Graph\n");
	start = high_resolution_clock::now();
	bool* graph = importGraph("ERgraphBig.txt", true, &n_nodes);
	stop = chrono::high_resolution_clock::now();
	printf("Done Importing Graph\n");
	duration = duration_cast<milliseconds>(stop - start);
	printf("Importing Graph: %.2d ms\n", duration);
	t2 = (int)max((double)t2, (double)n_nodes * 0.01);


	// Sequential BFS
	printf("\nStart Sequential BFS optimized\n");
	start = high_resolution_clock::now();
	lev = sequentialBasedBFS(graph, n_nodes, root);
	stop = chrono::high_resolution_clock::now();
	printf("Done Sequential BFS optimized\n");
	duration = duration_cast<milliseconds>(stop - start);
	printf("Sequential BFS optimized: %.2d ms\n", duration);

	plotLevelTable(lev, n_nodes);

	delete[] lev;


	// Queue Parallel BFS
	printf("\nStart Parallel BFS queue\n");
	start = high_resolution_clock::now();
	lev = queueBasedBFS(graph, n_nodes, root);
	stop = chrono::high_resolution_clock::now();
	printf("Done Parallel BFS queue\n");
	duration = duration_cast<milliseconds>(stop - start);
	printf("Parallel BFS queue: %.2d ms\n", duration);

	plotLevelTable(lev, n_nodes);

	delete[] lev;

	//printf("\nStart Parallel BFS read\n");
	//start = high_resolution_clock::now();
	//shared_lev = readBasedBFS();
	//stop = chrono::high_resolution_clock::now();
	//printf("Done Parallel BFS read\n");
	//duration = duration_cast<milliseconds>(stop - start);
	//printf("Parallel BFS read: %.2d us\n", duration);

	//plotLevelTable(shared_lev, shared_nodes);

	delete[] graph;
	delete[] lev;

	return 0;
}