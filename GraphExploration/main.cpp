#include "stdafx.h"

#include <iostream>
#include <iomanip>
#include <fstream>

bool* importGraph(const char* path, const bool undirected);
void printGraph(bool* graph);
int charsToInt(char chars[], int length);

int nodes = 0;

int main(int argc, char* argv[]) {

	bool* graph = importGraph("ERgraphTest.txt", true);

	printGraph(graph);

	return 0;
}

//int* levelSynchronousSequentialBFS(TUNGraph graph) {
//	const int nodes = graph.GetNodes();
//	int* output = new int[nodes];
//	TInt curr = 0;
//	output[0] = curr;
//	for (int i = 0; i < nodes; i++) {
//		curr = 
//	}
//}

bool* importGraph(const char* path, const bool undirected) {
	std::ifstream inFile;
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