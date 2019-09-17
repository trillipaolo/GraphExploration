#include "stdafx.h"
#include <string>

using namespace std;

void generateGraphs();

int main(int argc, char* argv[]) {
	generateGraphs();
	return 0;
}

void generateGraphs(){
	const int n_nodes = 93750;
	const int n_edges = 750000;

	const double a = 0.45;
	const double b = 0.25;
	const double c = 0.15;
 
	TInt::Rnd.PutSeed(0); // initialize random seed
	PUNGraph G;
	PNGraph Gd;
	int nodes;
	int edges;

	//nodes = n_nodes;
	//edges = n_edges;
	//for (int i = 0; i < 7; i++) {
	//	printf("Generating ER %d graph.\n", i);
	//	G = TSnap::GenRndGnm<PUNGraph>(nodes, edges, false);
	//	printf("Saving ER %d graph.\n", i);
	//	TSnap::SaveEdgeList(G, (string("ERcc") + to_string(i) + string(".txt")).c_str());
	//	printf("Done %d of 6.\n", i);

	//	nodes = floor(nodes / 2);
	//	edges = floor(edges / 4);
	//}

	//nodes = n_nodes;
	//edges = n_edges;
	//for (int i = 0; i < 7; i++) {
	//	printf("Generating RMAT %d graph.\n", i);
	//	Gd = TSnap::GenRMat(nodes, edges, a, b, c);
	//	printf("Saving RMAT graph.\n");
	//	TSnap::SaveEdgeList(Gd, (string("RMATcc") + to_string(i) + string(".txt")).c_str());
	//	printf("Done %d of 6.\n", i);

	//	nodes = floor(nodes / 2);
	//	edges = floor(edges / 4);
	//}


	//// Undirected Erdos - Renyi random graph.
	//printf("Generating ER small graph.\n");
	//G = TSnap::GenRndGnm<PUNGraph>(n_nodes/100, n_edges/100, false);
	//printf("Saving ER small graph.\n");
	//TSnap::SaveEdgeList(G, "ERgraphSmall.txt");
	//printf("Done 1 of 4.\n");

	printf("Generating ER big graph.\n");
	G = TSnap::GenRndGnm<PUNGraph>(n_nodes, n_edges, false);
	printf("Saving ER big graph.\n");
	TSnap::SaveEdgeList(G, "ERne05.txt");
	printf("Done 2 of 4.\n");

	//// Directed RMAT scale-free random graph.
	//printf("Generating RMAT small graph.\n");
	//Gd = TSnap::GenRMat(n_nodes/100, n_edges/100, a, b, c);
	//printf("Saving RMAT small graph.\n");
	//TSnap::SaveEdgeList(Gd, "RMATgraphSmall.txt");
	//printf("Done 3 of 4.\n");

	printf("Generating RMAT big graph.\n");
	Gd = TSnap::GenRMat(n_nodes, n_edges, a, b, c);
	printf("Saving RMAT big graph.\n");
	TSnap::SaveEdgeList(Gd, "RMATne05.txt");
	printf("Done 4 of 4.\n");
}
