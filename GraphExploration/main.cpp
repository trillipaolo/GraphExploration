#include <fstream>
#include <queue>
#include <chrono>
#include <thread>
#include <mutex>
#include <omp.h>


using namespace std;
using namespace std::chrono;

enum state { SEQ, QUEUE, QUEUE_TO_READ, READ, READ_TO_QUEUE, END };


// SHARED MEMORY //
int concurentThreadsSupported;
mutex mtx;

int t1 = 64;
int t2 = 262144;
int t3 = 2048;
float alpha = 2.0;
float beta = 2.0;
int queue_threshold = 20;


// FUNCTIONS IDENTIFIER //

short int* hybridBFS(bool* graph, const int n_nodes, const int root);
short int* sequentialBasedBFS(bool* graph, const int n_nodes, const int root);
short int* queueBasedBFS(bool* graph, const int n_nodes, const int root);
short int* readBasedBFS(bool* graph, int n_nodes, int root);

void plotLevelTable(short int* lev, int nodes);
queue<int> nextNodesQueue(bool* graph, int n_nodes, int curr);
vector<int> nextNodesVector(bool* graph, int n_nodes, int curr);
bool* importGraph(const char* path, const bool undirected, int* n_nodes);
bool* importGraphParallel(const char* path, const bool undirected, int* n_nodes);
void createGraph(bool* graph, int n_nodes, bool undirected, queue<int> &readArcs, bool &endOfFile);
void printGraph(bool* graph, int nodes);
int charsToInt(char chars[], int length);
bool scanBool();


// FUNCTIONS //

short int* hybridBFS(bool* graph, const int n_nodes, const int root) {
	omp_lock_t simple_lock;
	omp_init_lock(&simple_lock);

	state nextState = SEQ;
	bool exp = false;

	int n_curr = 0;
	int n_next = 0;

	vector<int> current;
	vector<int> next;

	// Initialize visited[] and lev[]
	short int level = 0;
	bool* visited = new bool[n_nodes];
	short int* lev = new short int[n_nodes];
	#pragma omp parallel for
	for (int i = 0; i < n_nodes; i++) {
		visited[i] = false;
		lev[i] = -1;
	}

	// initialize root node
	visited[root] = true;
	next.push_back(root);
	n_next++;
	lev[root] = 0;

	bool cond = true;
	while (cond) {
		cond = false;
		n_curr = n_next;
		n_next = 0;

		switch (nextState) {

		// SEQUENTIAL //
		case SEQ:

			current.clear();
			current.swap(next);

			for (int i = 0; i < current.size(); i++) {
				int curr = current.at(i);

				vector<int> neighbours = nextNodesVector(graph, n_nodes, curr);

				for (int j = 0; j < neighbours.size(); j++) {
					int neigh = neighbours.at(j);

					if (!visited[neigh]) {
						next.push_back(neigh);
						visited[neigh] = true;
						lev[neigh] = level + 1;
						cond = true;
					}
				}
			}

			n_next = next.size();

			if (n_next > t1) {
				nextState = QUEUE;
			}

			break;

		// QUEUE //
		case QUEUE:

			current.clear();
			current.swap(next);

			#pragma omp parallel
			{
				queue<int> private_queue;

				#pragma omp for
				for (int i = 0; i < current.size(); i++) {
					int curr = current.at(i);

					vector<int> neighbours = nextNodesVector(graph, n_nodes, curr);

					#pragma omp parallel for
					for (int j = 0; j < neighbours.size(); j++) {
						int neigh = neighbours.at(j);

						if (!visited[neigh]) {
							visited[neigh] = true;
							lev[neigh] = level + 1;
							cond = true;

							private_queue.push(neigh);

							if (private_queue.size() > queue_threshold) {
								if (omp_test_lock(&simple_lock)) {
									while (!private_queue.empty()) {
										next.push_back(private_queue.front());
										private_queue.pop();
									}
									omp_unset_lock(&simple_lock);
								}
							}
						}
					}
				}
				omp_set_lock(&simple_lock);
				while (!private_queue.empty()) {
					next.push_back(private_queue.front());
					private_queue.pop();
				}
				omp_unset_lock(&simple_lock);
			}

			n_next = next.size();

			if (n_next > t2 || (n_next > alpha * n_curr)) {
				nextState = QUEUE_TO_READ;
				exp = n_next > alpha * n_curr;
			}
			else if (n_next < t1) {
				nextState = SEQ;
			}
			break;

		// QUEUE to READ
		case QUEUE_TO_READ:

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
						cond = true;

						#pragma omp critical (nextPush)
						{
							n_next++;
						}
					}
				}
			}

			nextState = READ;
			break;

		// READ //
		case READ:

			#pragma omp parallel for
			for (int curr = 0; curr < n_nodes; curr++) {
				if (lev[curr] == level) {

					vector<int> neighbours = nextNodesVector(graph, n_nodes, curr);

					#pragma omp parallel for
					for (int j = 0; j < neighbours.size(); j++) {
						int neigh = neighbours.at(j);

						if (!visited[neigh]) {
							visited[neigh] = true;
							lev[neigh] = level + 1;
							cond = true;

							#pragma omp critical (nextPush)
							{
								n_next++;
							}
						}
					}
				}
			}

			if (n_next < t2 || (exp && (n_next < beta * n_curr))) {
				nextState = READ_TO_QUEUE;
			}
			break;

		// READ to QUEUE //
		case READ_TO_QUEUE:

			next.clear();

			#pragma omp parallel
			{
				queue<int> private_queue;

				#pragma omp for
				for (int curr = 0; curr < n_nodes; curr++) {
					if (lev[curr] == level) {

						vector<int> neighbours = nextNodesVector(graph, n_nodes, curr);

						#pragma omp parallel for
						for (int j = 0; j < neighbours.size(); j++) {
							int neigh = neighbours.at(j);

							if (!visited[neigh]) {
								visited[neigh] = true;
								lev[neigh] = level + 1;
								cond = true;

								private_queue.push(neigh);

								if (private_queue.size() > queue_threshold) {
									if (omp_test_lock(&simple_lock)) {
										while (!private_queue.empty()) {
											next.push_back(private_queue.front());
											private_queue.pop();
										}
										omp_unset_lock(&simple_lock);
									}
								}
							}
						}
					}
				}
				omp_set_lock(&simple_lock);
				while (!private_queue.empty()) {
					next.push_back(private_queue.front());
					private_queue.pop();
				}
				omp_unset_lock(&simple_lock);
			}

			n_next = next.size();

			if (n_next > t1) {
				nextState = QUEUE;
			}
			else {
				nextState = SEQ;
			}

			break;
		default:
			break;
		}

		level++;
	}

	delete[] visited;

	omp_destroy_lock(&simple_lock);

	return lev;
}

short int* sequentialBasedBFS(bool* graph, const int n_nodes, const int root) {

	queue<int> current;
	queue<int> next;
	queue<int> neighbours;

	// Initialize visited[] and lev[]
	short int level = 0;
	bool* visited = new bool[n_nodes];
	short int* lev = new short int[n_nodes];
	for (int i = 0; i < n_nodes; i++) {
		visited[i] = false;
		lev[i] = -1;
	}

	// Initialize root node
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

short int* queueBasedBFS(bool* graph, const int n_nodes, const int root) {
	omp_lock_t simple_lock;
	omp_init_lock(&simple_lock);

	vector<int> current;
	vector<int> next;

	// Initialize visited[] and lev[]
	short int level = 0;
	bool* visited = new bool[n_nodes];
	short int* lev = new short int[n_nodes];
	#pragma omp parallel for
	for (int i = 0; i < n_nodes; i++) {
		visited[i] = false;
		lev[i] = -1;
	}

	// initialize root node
	visited[root] = true;
	next.push_back(root);
	lev[root] = 0;

	while (!next.empty()) {

		current.clear();
		current.swap(next);

		#pragma omp parallel
		{
			queue<int> private_next;

			#pragma omp for
			for (int i = 0; i < current.size(); i++) {
				int curr = current.at(i);

				vector<int> neighbours = nextNodesVector(graph, n_nodes, curr);

				#pragma omp parallel for
				for (int j = 0; j < neighbours.size(); j++) {
					int neigh = neighbours.at(j);
					if (!visited[neigh]) {
						visited[neigh] = true;
						lev[neigh] = level + 1;

						private_next.push(neigh);

						if (private_next.size() > queue_threshold) {
							if (omp_test_lock(&simple_lock)) {
								while (!private_next.empty()) {
									next.push_back(private_next.front());
									private_next.pop();
								}
								omp_unset_lock(&simple_lock);
							}
						}
					}
				}
			}

			omp_set_lock(&simple_lock);
			while (!private_next.empty()) {
				next.push_back(private_next.front());
				private_next.pop();
			}
			omp_unset_lock(&simple_lock);
		}

		level++;
	}

	delete[] visited;

	omp_destroy_lock(&simple_lock);

	return lev;
}

short int* readBasedBFS(bool* graph, int n_nodes, int root) {
	
	// Initialize visited[] and lev[]
	short int level = 0;
	bool* visited = new bool[n_nodes];
	short int* lev = new short int[n_nodes];
	#pragma omp parallel for
	for (int i = 0; i < n_nodes; i++) {
		visited[i] = false;
		lev[i] = -1;
	}

	// Initialize root node
	visited[root] = true;
	lev[root] = 0;

	bool cond = true;
	while (cond) {
		cond = false;

		#pragma omp parallel for
		for (int curr = 0; curr < n_nodes; curr++) {
			if (lev[curr] == level) {
				
				vector<int> neighbours = nextNodesVector(graph, n_nodes, curr);

				#pragma omp parallel for
				for (int j = 0; j < neighbours.size(); j++) {
					int neigh = neighbours.at(j);

					if (!visited[neigh]) {
						visited[neigh] = true;
						lev[neigh] = level + 1;
						cond = true;
					}
				}
			}
		}

		level++;
	}

	delete[] visited;

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

	printf("___________________________\n");
	printf("\n Level\tNodes\tPercentage\n");
	printf("___________________________\n\n");
	for (int i = 1; i < maxLev + 2; i++) {
		printf(" %d\t%d\t%.2f %%\n", i - 1, count[i], (float)((count[i]) * 100) / nodes);
	}
	printf("___________________________\n\n");
	printf(" Total:\t%d\t%.2f %%\n", nodes - count[0], (float)((nodes - count[0]) * 100) / nodes);
	printf("___________________________\n\n");

	delete[] count;
}

queue<int> nextNodesQueue(bool* graph, int n_nodes, int curr) {

	long long int lc = (long long int)curr;
	long long int ln_ns = (long long int)n_nodes;
	queue<int> next;
	for (long long int i = 0; i < ln_ns; i++) {
		if (graph[lc * ln_ns + i]) {
			next.push((int)i);
		}
	}

	return next;
}

vector<int> nextNodesVector(bool* graph, int n_nodes, int curr)  {

	long long int lc = (long long int)curr;
	long long int ln_ns = (long long int)n_nodes;
	vector<int> neighbours;


	// WARNING: USE ONLY ONE METHOD BETWEEN LOCK AND PRIVATE QUEUE, COMMENT THE OTHER ONE

	//// LOCK:
	////	lighter version, in my case works better
	omp_lock_t neigh_lock;
	omp_init_lock(&neigh_lock);
	#pragma omp parallel for
	for (long long int i = 0; i < ln_ns; i++) {
		if (graph[lc * ln_ns + i]) {
			omp_set_lock(&neigh_lock);
			neighbours.push_back((int)i);
			omp_unset_lock(&neigh_lock);
		}
	}
	omp_destroy_lock(&neigh_lock);


	//// PRIVATE QUEUE:
	////		to not stop execution of ohter threads at every vector push_back
	//omp_lock_t neigh_lock;
	//omp_init_lock(&neigh_lock);
	//#pragma omp parallel
	//{
	//	queue<int> private_neighbours;
	//	#pragma omp for
	//	for (long long int i = 0; i < ln_ns; i++) {
	//		if (graph[lc * ln_ns + i]) {
	//			private_neighbours.push((int)i);
	//		}
	//		if (private_neighbours.size() > queue_threshold) {
	//			if (omp_test_lock(&neigh_lock)) {
	//				while (!private_neighbours.empty()) {
	//					neighbours.push_back(private_neighbours.front());
	//					private_neighbours.pop();
	//				}
	//				omp_unset_lock(&neigh_lock);
	//			}
	//		}
	//	}
	//	omp_set_lock(&neigh_lock);
	//	while (!private_neighbours.empty()) {
	//		neighbours.push_back(private_neighbours.front());
	//		private_neighbours.pop();
	//	}
	//	omp_unset_lock(&neigh_lock);
	//}
	//omp_destroy_lock(&neigh_lock);


	return neighbours;
}

bool* importGraph(const char* path, const bool undirected, int* n_nodes) {
	ifstream inFile;
	inFile.open(path);
	if (!inFile) {
		throw "File not Found!";
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

bool* importGraphParallel(const char* path, const bool undirected, int* n_nodes) {
	ifstream inFile;
	inFile.open(path);
	if (!inFile) {
		throw "File not Found!";
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
	long long int ln = (long long int) * n_nodes;

	// Arrive at edges
	inFile.getline(curr, 256, '#');
	inFile.getline(curr, 256, '\n');

	// Initialize graph without edges
	bool* graph;
	try {
		graph = new bool[ln * ln];
	}
	catch (const std::bad_alloc& e) {
		printf("Allocation failed: %s of bool[%lld]\n", e.what(), ln*ln);
		inFile.close();
		exit(-1);
	}

	#pragma omp parallel for
	for (long long int i = 0; i < ln * ln; i++) {
		graph[i] = false;
	}

	// Read edges
	bool endOfFile = false;
	queue<int> readArcs;
	int sn;
	int en;

	thread* at = new thread[concurentThreadsSupported - 1];
	for (int i = 0; i < (concurentThreadsSupported - 1); i++) {
		at[i] = thread(createGraph, graph, *n_nodes, undirected, ref(readArcs), ref(endOfFile));
	}

	while (!inFile.getline(curr, 256, '\t').eof()) {
		sn = charsToInt(curr, 256);

		inFile.getline(curr, 256, '\n');
		en = charsToInt(curr, 256);

		mtx.lock();
		readArcs.push(sn);
		readArcs.push(en);
		mtx.unlock();
	}

	endOfFile = true;

	inFile.close();

	for (int i = 0; i < (concurentThreadsSupported - 1); i++) {
		at[i].join();
	}

	delete[] at;

	return graph;
}

void createGraph(bool* graph, int n_nodes, bool undirected, queue<int> &readArcs, bool &endOfFile) {
	long long int ln = n_nodes;
	long long int sn;
	long long int en;
	bool terminate = false;

	while (!terminate) {
		mtx.lock();
		if (!readArcs.empty()) {
			sn = readArcs.front();
			readArcs.pop();
			en = readArcs.front();
			readArcs.pop();
			mtx.unlock();

			graph[sn * ln + en] = true;
			if (undirected) {
				graph[en * ln + sn] = true;
			}
		}
		else {
			mtx.unlock();
			if (endOfFile) {
				terminate = true;
			}
			else {
				this_thread::sleep_for(milliseconds(100));
			}
		}
	}
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

bool scanBool() {
	char input;
	do {
		scanf_s("%c", &input, 1);
	} while (input == '\n');
	return input == '1';
}


// MAIN //

int main(int argc, char* argv[]) {
	concurentThreadsSupported = thread::hardware_concurrency();

	printf("Decide which algorthm to run:\n");
	printf(" - S: Sequential\n");
	printf(" - Q: Queue\n");
	printf(" - R: Read\n");
	printf(" - H: Hybrid\n");
	char algorithm;
	do {
		printf("Insert a character: ");
		scanf_s("%c", &algorithm, 1);
		if (algorithm >= 'a' && algorithm <= 'z') {
			algorithm = algorithm - 'a' + 'A';
		}
	} while (algorithm != 'S' && algorithm != 'Q' && algorithm != 'R' && algorithm != 'H');

	printf("\n");
	int times;
	do {
		string input;
		printf("Decide how many runs : ");
		scanf_s("%s", &input, 10);
		times = atoi(input.c_str());
		if (times <= 0) {
			printf("Insert a positive number!\n");
		}
		else if (times > 50) {
			printf("If you run it more than 50 times my laptop will explode :(\n");
		}
	} while (times <= 0 || times > 50);


	// For time computation
	time_point<steady_clock> start;
	time_point<steady_clock> stop;
	milliseconds* duration = new milliseconds[10];


	// Import graph
	bool pathFound;
	string path;
	const int root = 0;
	short int* lev;
	int n_nodes = -1;
	bool* graph;
	do {
		pathFound = true;
		printf("\nInsert the name of file's graph: ");
		scanf_s("%s", &path, 50);
		printf("The graph is undirected (0|1)? ");
		bool undirected = scanBool();
		printf("\nStart Importing Graph\n");
		start = high_resolution_clock::now();
		try {
			graph = importGraphParallel(path.c_str(), undirected, &n_nodes);
		}
		catch (const char* msg) {
			pathFound = false;
			graph = NULL;
			printf(msg);
			printf("\n");
		}
		stop = chrono::high_resolution_clock::now();
		printf("Done Importing Graph\n");
		duration[0] = duration_cast<milliseconds>(stop - start);
		printf("Importing Graph: %d ms\n\n", (int)duration[0].count());
		t2 = (int)max((double)t2, (double)n_nodes * 0.01);
	} while (!pathFound);


	// Run choosen algorithm
	for (int i = 0; i < times; i++) {

		printf("Start BFS: %d of %d\n", i+1, times);
		start = high_resolution_clock::now();
		switch (algorithm) {
		case 'H':
			lev = hybridBFS(graph, n_nodes, root);
			break;
		case 'R':
			lev = readBasedBFS(graph, n_nodes, root);
			break;
		case 'Q':
			lev = queueBasedBFS(graph, n_nodes, root);
			break;
		case 'S':
			lev = sequentialBasedBFS(graph, n_nodes, root);
			break;
		default:
			lev = NULL;
			break;
		}
		stop = chrono::high_resolution_clock::now();
		duration[i] = duration_cast<milliseconds>(stop - start);
		printf("Duration: %d ms\n", (int)duration[i].count());

		if (i == times - 1) {
			printf("\n");
			plotLevelTable(lev, n_nodes);
		}

		delete[] lev;
	}


	// Print results
	printf("\n");
	switch (algorithm) {
	case 'H':
		printf("Hybrid");
		break;
	case 'R':
		printf("Read");
		break;
	case 'Q':
		printf("Queue");
		break;
	case 'S':
		printf("Sequential");
		break;
	default:
		break;
	}
	printf(" BFS:\nIter\tTime (ms)\n");
	int mean = 0;
	for (int i = 0; i < times; i++) {
		mean = mean + (int)duration[i].count();
		printf("%d\t%d\n", i, (int)duration[i].count());
	}
	mean = mean / times;
	printf("___________________\n\nMean\t%d\n", mean);


	delete[] graph;


	// Needed to be able to read the results
	string noclose;
	scanf_s("%c", noclose, 50);
	scanf_s("%c", noclose, 50);

	return 0;
}