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

int t1 = 64;
int t2 = 262144;
int t3 = 2048;
float alpha = 2.0;
float beta = 2.0;

int shared_nodes = 0;
bool* shared_graph;
short int* shared_lev;
bool shared_cond = false;
int shared_level = 0;


// FUNCTIONS IDENTIFIER

short int* levelSynchronousSequentialBFSbasic();

short int* sequentialBasedBFS();

short int* readBasedBFS();
void readBFSarray(int num);
void readBFScurrent(int curr);

short int* queueBasedBFS();
void queueBFScurrent(int curr, bool* visited, queue<int>& next);
void queueBFSneighbours(int neigh, bool* visited, queue<int>& next);

void initializeVisitedLev(int num, bool* visited, short int* lev);
void plotLevelTable(short int* lev, int nodes);
queue<int> nextNodesQueue(bool* graph, int node, int nodes);
set<int> nextNodesSet(bool* graph, int node, int nodes);
bool* importGraph(const char* path, const bool undirected);
void printGraph(bool* graph, int nodes);
int charsToInt(char chars[], int length);


// FUNCTIONS

short int* hybridBFS() {
	state nextState = SEQ;

	int n_curr;
	int n_next;
	bool exp = false;

	shared_cond = true;
	while (shared_cond) {
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
}

short int* levelSynchronousSequentialBFSbasic() {
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
	short int* lev = new short int[shared_nodes];
	for (int i = 0; i < shared_nodes; i++) {
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

			set<int> neighbours = nextNodesSet(shared_graph, *itr, shared_nodes);
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

short int* sequentialBasedBFS() {
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

	shared_level = 0;
	bool* visited = new bool[shared_nodes];

	thread* at = new thread[concurentThreadsSupported];
	for (int i = 0; i < concurentThreadsSupported; i++) {
		at[i] = thread(initializeVisitedLev, i, visited, shared_lev);
	}
	for (int i = 0; i < concurentThreadsSupported; i++) {
		at[i].join();
	}

	queue<int> current;
	queue<int> next;
	queue<int> neighbours;

	// root is node 0
	visited[0] = true;
	next.push(0);
	shared_lev[0] = 0;

	int curr;
	int neigh;
	while (!next.empty()) {
		// current = queue<int>();		clean a queue
		swap(current,next);

		while (!current.empty()) {
			curr = current.front();

			neighbours = nextNodesQueue(shared_graph, curr, shared_nodes);
			while (!neighbours.empty()) {
				neigh = neighbours.front();

				if (!visited[neigh]) {
					next.push(neigh);
					visited[neigh] = true;
					shared_lev[neigh] = shared_level + 1;
				}

				neighbours.pop();
			}

			current.pop();
		}
		shared_level++;
	}

	delete[] at;
	delete[] visited;

	return shared_lev;
}

short int* readBasedBFS() {
	//BFS_Read(G: Graph, r : Node) {
	//	Bitmap V;
	//	Bool fin[threads];
	//	V.set(r.id);
	//	int level = 0; r.lev = level;
	//	bool finished = false;
	//	while (!finished) {
	//		fork;
	//		fin[tid] = true;
	//		foreach(c: G.Nodes.partition(tid)) {
	//			if (c.lev != level)		// if it is not current, go to the next node
	//				continue;
	//			foreach(n: c.nbrs) {
	//				if (!V.isSet(n.id)) { // test and test-and-set
	//					if (V.atomicSet(n.id)) {
	//						n.lev = level + 1;
	//						fin[tid] = false;
	//					}
	//				}
	//			}
	//		}
	//		join;
	//		finished = logicalAnd(fin, threads);
	//		level++;
	//	}
	//}

	queue<thread> threads;
	
	shared_level = 0;
	bool* visited = new bool[shared_nodes];

	thread* at = new thread[concurentThreadsSupported];
	for (int i = 0; i < concurentThreadsSupported; i++) {
		at[i] = thread(initializeVisitedLev, i, visited, shared_lev);
	}
	for (int i = 0; i < concurentThreadsSupported; i++) {
		at[i].join();
	}

	// root is node 0
	visited[0] = true;
	shared_lev[0] = 0;

	shared_cond = true;
	while (shared_cond) {
		shared_cond = false;
		for (int i = 0; i < concurentThreadsSupported; i++) {
			at[i] = thread(readBFSarray, i);
		}
		for (int i = 0; i < concurentThreadsSupported; i++) {
			at[i].join();
		}

		//for (int i = 0; i < nodes; i++) {
		//	if (lev[i] == level) {
		//		cond = true;
		//		threads.push(thread(readBFScurrent, graph, nodes, i, lev, level));
		//	}
		//}

		//while (!threads.empty()) {
		//	threads.front().join();
		//	threads.pop();
		//}

		shared_level++;
	}

	delete[] at;
	delete[] visited;

	return shared_lev;
}

void readBFSarray(int num) {
	int start = (shared_nodes * num) / concurentThreadsSupported;
	int end = (shared_nodes * (num + 1)) / concurentThreadsSupported;
	if (end > shared_nodes) {
		end = shared_nodes;
	}

	queue<thread> threads;

	for (int i = start; i < end; i++) {
		if (shared_lev[i] == shared_level) {
			shared_cond = true;
			threads.push(thread(readBFScurrent, i));
		}
	}
		
	while (!threads.empty()) {
		threads.front().join();
		threads.pop();
	}
}

void readBFScurrent(int curr) {
	queue<int> neighbours;
	neighbours = nextNodesQueue(shared_graph, curr, shared_nodes);

	int neigh;
	while (!neighbours.empty()) {

		neigh = neighbours.front();

		if (shared_lev[neigh] == -1) {
			shared_lev[neigh] = shared_level + 1;
		}

		neighbours.pop();
	}
}

short int* queueBasedBFS() {

	queue<int> current;
	queue<int> next;
	queue<thread> threads;

	shared_level = 0;
	bool* visited = new bool[shared_nodes];

	thread* at = new thread[concurentThreadsSupported];
	for (int i = 0; i < concurentThreadsSupported; i++) {
		at[i] = thread(initializeVisitedLev, i, visited, shared_lev);
	}
	for (int i = 0; i < concurentThreadsSupported; i++) {
		at[i].join();
	}

	// root is node 0
	visited[0] = true;
	next.push(0);
	shared_lev[0] = 0;

	int curr;
	while (!next.empty()) {
		swap(current, next);

		while (!current.empty()) {

			curr = current.front();

			threads.push(thread(queueBFScurrent, curr, visited, ref(next)));

			current.pop();
		}

		while (!threads.empty()) {
			threads.front().join();
			threads.pop();
		}

		shared_level++;
	}

	delete[] at;
	delete[] visited;

	return shared_lev;
}

void queueBFScurrent(int curr, bool* visited, queue<int> &next) {
	queue<int> neighbours;
	queue<thread> threads;

	int neigh;

	neighbours = nextNodesQueue(shared_graph, curr, shared_nodes);
	while (!neighbours.empty()) {

		neigh = neighbours.front();

		if (!visited[neigh]) {
			threads.push(thread(queueBFSneighbours, neigh, visited, ref(next)));
		}

		neighbours.pop();
	}

	while (!threads.empty()) {
		threads.front().join();
		threads.pop();
	}
}

void queueBFSneighbours(int neigh, bool* visited, queue<int> &next) {
	if (!visited[neigh]) {
		visited[neigh] = true;
		shared_lev[neigh] = shared_level + 1;

		mtx.lock();
		next.push(neigh);
		mtx.unlock();
	}
}

void initializeVisitedLev(int num, bool* visited, short int* lev) {
	int start = (shared_nodes * num) / concurentThreadsSupported;
	int end = (shared_nodes * (num + 1)) / concurentThreadsSupported;
	if (end > shared_nodes) {
		end = shared_nodes;
	}
	for (int i = start; i < end; i++) {
		visited[i] = false;
		lev[i] = -1;
	}
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

	delete[] count;
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
	shared_nodes = charsToInt(curr, 256);
	t2 = (int) max((double)t2, (double)shared_nodes * 0.01);
	long long int ln = (long long int)shared_nodes;

	// Arrive at edges
	inFile.getline(curr, 256, '#');
	inFile.getline(curr, 256, '\n');

	// Initialize graph without edges
	bool* graph = new bool[ln * ln];
	for (int i = 0; i < shared_nodes * shared_nodes; i++) {
		graph[i] = false;
	}

	// Read edges
	long long int sn;
	long long int en;
	while (!inFile.getline(curr, 256, '\t').eof()) {
		sn = charsToInt(curr, 256);

		inFile.getline(curr, 256, '\n');
		en = charsToInt(curr, 256);

		graph[sn * shared_nodes + en] = true;
		if (undirected) {
			graph[en * shared_nodes + sn] = true;
		}
	}

	inFile.close();

	return graph;
}

void printGraph(bool* graph, int nodes) {
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


// MAIN

int main(int argc, char* argv[]) {

	omp_set_dynamic(0);

	#pragma omp parallel num_threads(5)
		{
			printf("Hello World... from thread = %d\n", omp_get_thread_num());
		}

	concurentThreadsSupported = thread::hardware_concurrency();

	printf("MAX_THREADS: %d\n", concurentThreadsSupported);

	time_point<steady_clock> start;
	time_point<steady_clock> stop;
	milliseconds duration;

	shared_graph = importGraph("ERgraphBig.txt", true);
	shared_lev = new short int[shared_nodes];

	////printGraph(shared_graph, shared_nodes);

	////printf("\nStart Sequential BFS normal\n");
	////start = high_resolution_clock::now();
	////shared_lev = levelSynchronousSequentialBFSbasic();
	////stop = chrono::high_resolution_clock::now();
	////printf("Done Sequential BFS normal\n");
	////duration = duration_cast<milliseconds>(stop - start);
	////printf("Sequential BFS basic: %.2d us\n", duration);

	////plotLevelTable(shared_lev, shared_nodes);

	////printf("\nStart Sequential BFS optimized\n");
	////start = high_resolution_clock::now();
	////shared_lev = sequentialBasedBFS();
	////stop = chrono::high_resolution_clock::now();
	////printf("Done Sequential BFS optimized\n");
	////duration = duration_cast<milliseconds>(stop - start);
	////printf("Sequential BFS optimized: %.2d us\n", duration);

	////plotLevelTable(shared_lev, shared_nodes);


	//printf("\nStart Parallel BFS read\n");
	//start = high_resolution_clock::now();
	//shared_lev = readBasedBFS();
	//stop = chrono::high_resolution_clock::now();
	//printf("Done Parallel BFS read\n");
	//duration = duration_cast<milliseconds>(stop - start);
	//printf("Parallel BFS read: %.2d us\n", duration);

	//plotLevelTable(shared_lev, shared_nodes);


	//printf("\nStart Parallel BFS queue\n");
	//start = high_resolution_clock::now();
	//shared_lev = queueBasedBFS();
	//stop = chrono::high_resolution_clock::now();
	//printf("Done Parallel BFS queue\n");
	//duration = duration_cast<milliseconds>(stop - start);
	//printf("Parallel BFS queue: %.2d us\n", duration);

	//plotLevelTable(shared_lev, shared_nodes);

	//delete[] shared_graph;
	//delete[] shared_lev;

	return 0;
}