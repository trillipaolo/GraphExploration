#include "stdafx.h"

void generateGraphs();

int main(int argc, char* argv[]) {
	generateGraphs();
	return 0;
}

void generateGraphs(){
	const int n_nodes = 120000;
	const int n_edges = 2560000;

	const double a = 0.45;
	const double b = 0.25;
	const double c = 0.15;
 
	TInt::Rnd.PutSeed(0); // initialize random seed
	PUNGraph G;
	PNGraph Gd;

	//// Undirected Erdos - Renyi random graph.
	//printf("Generating ER small graph.\n");
	//G = TSnap::GenRndGnm<PUNGraph>(n_nodes/100, n_edges/100, false);
	//printf("Saving ER small graph.\n");
	//TSnap::SaveEdgeList(G, "ERgraphSmall.txt");
	//printf("Done 1 of 4.\n");

	//printf("Generating ER big graph.\n");
	//G = TSnap::GenRndGnm<PUNGraph>(n_nodes, n_edges, false);
	//printf("Saving ER big graph.\n");
	//TSnap::SaveEdgeList(G, "ERgraphBig.txt");
	//printf("Done 2 of 4.\n");

	//// Directed RMAT scale-free random graph.
	//printf("Generating RMAT small graph.\n");
	//Gd = TSnap::GenRMat(n_nodes/100, n_edges/100, a, b, c);
	//printf("Saving RMAT small graph.\n");
	//TSnap::SaveEdgeList(Gd, "RMATgraphSmall.txt");
	//printf("Done 3 of 4.\n");

	//printf("Generating RMAT big graph.\n");
	//Gd = TSnap::GenRMat(n_nodes, n_edges, a, b, c);
	//printf("Saving RMAT big graph.\n");
	//TSnap::SaveEdgeList(Gd, "RMATgraphBig.txt");
	//printf("Done 4 of 4.\n");

	printf("Generating ER big graph.\n");
	G = TSnap::GenRndGnm<PUNGraph>(n_nodes / 10, n_edges / 10, false);
	printf("Saving ER big graph.\n");
	TSnap::SaveEdgeList(G, "ERgraphMedium.txt");
	printf("Done 2 of 4.\n");

	// Directed RMAT scale-free random graph.
	printf("Generating RMAT small graph.\n");
	Gd = TSnap::GenRMat(n_nodes / 10, n_edges / 10, a, b, c);
	printf("Saving RMAT small graph.\n");
	TSnap::SaveEdgeList(Gd, "RMATgraphMedium.txt");
	printf("Done 3 of 4.\n");
}
