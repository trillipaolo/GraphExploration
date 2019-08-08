#include "stdafx.h""

int main(int argc, char* argv[]) {
	PUNGraph G = TSnap::LoadEdgeList<PUNGraph>("ERgraphSmall.txt", 0, 1);
	printf("Hello World");
	return 0;
}