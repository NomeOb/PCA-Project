// Copyright 2016 Mohamed Essam
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <utility>
#include <vector>
#include <map>
#include "lodepng.h"

#define PI 3.141569

struct Edge{
	unsigned int from, to, weight;
};

struct Graph{
	Edge* edges;
	int edgeCount;
};

// Graph Creation
int square(int a) {
	return a*a;
}

int getEuclidianDistance(int r1, int g1, int b1, int r2, int g2, int b2) {
	return 4 * sqrt((double)(square(r1 - r2) + square(g1 - g2) + square(b1 - b2)));
}

Graph* createGraph(unsigned char* rgb, unsigned int maxWeight,
	unsigned int height, unsigned int width) {
	unsigned int edgeCount = 2 * width * height - width - height;
	Edge* graph = new Edge[edgeCount];
	unsigned int nextEdge = 0;
	unsigned int histogram[maxWeight];
	memset(histogram, 0, sizeof(histogram));
	edgeCount = 0;
	int a = 0, b = 0;
	for (int j = 0; j < width; j++) {
		for (int i = 0; i < height; i++) {
			unsigned int idx1, idx2, idx3;
			idx1 = (i*width + j) * 3;
			idx3 = ((i + 1)*width + j) * 3;
			if (i + 1 < height) {
				unsigned int w = getEuclidianDistance(rgb[idx1], rgb[idx1 + 1],
					rgb[idx1 + 2], rgb[idx3], rgb[idx3 + 1], rgb[idx3 + 2]);
				//printf("n1 %d n2 %d 2 %d\n", idx1/3, idx3/3, w);
				if (w <= maxWeight) {
					graph[nextEdge++] = { idx1 / 3, idx3 / 3, w };
					histogram[w]++;
					edgeCount++;
				}
				else a++;
			}
			idx2 = (i*width + j + 1) * 3;
			if (j + 1 < width) {
				unsigned int w = getEuclidianDistance(rgb[idx1], rgb[idx1 + 1],
					rgb[idx1 + 2], rgb[idx2], rgb[idx2 + 1], rgb[idx2 + 2]);
				//printf("n1 %d n2 %d 2 %d\n", idx1/3, idx2/3, w);
				if (w <= maxWeight) {
					graph[nextEdge++] = { idx1 / 3, idx2 / 3, w };
					histogram[w]++;
					edgeCount++;
				}
				else b++;
			}
		}
	}
	printf("%d %d", a, b);
	for (int i = 1; i < maxWeight; i++) {
		histogram[i] += histogram[i - 1];
	}
	Edge* sortedGraph = new Edge[edgeCount];
	int st_pos = 0;
	for (int i = 0; i < edgeCount; i++) {
		int w = graph[i].weight;
		if (w == 0) {
			sortedGraph[st_pos++] = graph[i];
		}
		else {
			sortedGraph[histogram[w - 1]++] = graph[i];
		}
	}
	delete[] graph;
	Graph* ret = new Graph();
	ret->edgeCount = edgeCount;
	ret->edges = sortedGraph;
	return ret;
}

// Initial Regions
struct UnionFind{
	unsigned int* parent;
	unsigned int* rank;

	explicit UnionFind(unsigned int size) {
		parent = new unsigned int[size];
		rank = new unsigned int[size];
		for (unsigned int i = 0; i < size; i++) {
			rank[i] = 1;
			parent[i] = i;
		}
	}

	~UnionFind() {
		delete[] parent;
		delete[] rank;
	}

	unsigned int find(unsigned int a) {
		return parent[a] = (a == parent[a] ? a : find(parent[a]));
	}

	unsigned int unite(unsigned int a, unsigned int b) {
		a = find(a), b = find(b);
		if (a == b) return a;
		if (rank[b] > rank[a]) std::swap(a, b);
		parent[b] = a;
		rank[a] += rank[b];
		return a;
	}

	bool isRepresetative(unsigned int x) {
		return x == find(x);
	}

	bool isUnited(unsigned int a, unsigned int b) {
		return find(a) == find(b);
	}
};

struct Region{
	unsigned int representative, size, credit;
};

struct RegionalData {
	UnionFind* uf;
	Region* regions;

	RegionalData(UnionFind* _uf, Region* _regions) {
		uf = _uf;
		regions = _regions;
	}
};

int computeCredit(int regionSize) {
	return sqrt(4 * PI * regionSize);
}

RegionalData* getInitialRegions(Edge* graph, int minWeight, int minRegionSize,
	int edgeCount, int width, int height) {
	UnionFind* uf = new UnionFind(width*height);
	for (int i = 0; i < edgeCount; i++) {
		if (graph[i].weight <= minWeight){
			uf->unite(graph[i].from, graph[i].to);
		}
	}
	// for (int i = 0; i < edgeCount && graph[i].weight <= minWeight; i++) {
	// 	if (graph[i].weight > minWeight)
	// 		continue;
	// 	uf->unite(graph[i].from, graph[i].to);
	// }
	Region* regions = new Region[width*height];
	for (int i = 0; i < width*height; i++) {
		regions[i].size = regions[i].representative = regions[i].credit = 0;
	}
	for (int i = 0; i < width*height; i++) {
		regions[uf->find(i)].size++;
		regions[uf->find(i)].representative = uf->find(i);
	}
	for (int i = 0; i < width*height; i++) {
		if (regions[i].size >= minRegionSize)
			regions[i].credit = computeCredit(regions[i].size);
		else{
			regions[i].credit = 100000;
		}
	}
	RegionalData* ret = new RegionalData(uf, regions);
	return ret;
}

// Region Expansion

struct SegmentationData{
	int* segment;
	int segmentCount;
};

SegmentationData expandRegions(RegionalData* rd, int minWeight,
	int width, int height, Edge* graph, int edgeCount) {
	int united_cnt = 0;
	for (int i = 0; i < edgeCount; i++) {
		if (graph[i].weight <= minWeight) continue;
		if (!rd->uf->isUnited(graph[i].from, graph[i].to)) {
			int a = rd->uf->find(graph[i].from);
			int b = rd->uf->find(graph[i].to);
			int credit = std::min(rd->regions[a].credit, rd->regions[b].credit);
			if (credit > graph[i].weight) {
				int s = rd->uf->unite(a, b);
				rd->regions[s].credit = credit - graph[i].weight;
				united_cnt++;
			}
		}
	}
	printf("United cnt%d\n", united_cnt);
	int* repTable = new int[width*height];
	int idx = 0;
	for (int i = 0; i < width*height; i++) {
		if (rd->uf->isRepresetative(i)) {
			repTable[i] = idx++;
		}
	}
	int* regionals = new int[width*height];
	for (int i = 0; i < width*height; i++) {
		regionals[i] = repTable[rd->uf->find(i)];
	}
	SegmentationData ret = { regionals, idx };
	return ret;
}

// Main entry point

SegmentationData segmentImage(unsigned char* rgb, int height,
	int width, int minWeight, int maxWeight, int minRegionSize) {
	Graph* graph = createGraph(rgb, maxWeight, height, width);
	RegionalData* regional = getInitialRegions(graph->edges, minWeight,
		minRegionSize, graph->edgeCount, width, height);
	return expandRegions(regional, minWeight, width, height,
		graph->edges, graph->edgeCount);
}

int main(){
	std::vector<unsigned char> tmp;
	unsigned int w, h;
	lodepng::decode(tmp, w, h, "4squaresSmall.png");
	unsigned char* rgb = new unsigned char[3 * w*h];
	for (int i = 0; i < w*h; i++)
	{
		rgb[i * 3] = tmp[i * 4];
		rgb[i * 3 + 1] = tmp[i * 4 + 1];
		rgb[i * 3 + 2] = tmp[i * 4 + 2];
		//printf("\n");
		//printf("Pixel at [%d]: R = %u G= %u B= %u \n",i, rgb[i * 3], rgb[i * 3 + 1] , rgb[i * 3 +2]);
	} 
	SegmentationData data = segmentImage(rgb, h, w, 20, 90, 50);
	std::vector<int> v(data.segmentCount, 0);
	for (int i = 0; i < w*h; i++)
	{
		v[data.segment[i]]++;
	}
	std::vector<int> segColorsR(data.segmentCount, 0);
	std::vector<int> segColorsG(data.segmentCount, 0);
	std::vector<int> segColorsB(data.segmentCount, 0);
	int idx = 0;
	std::vector<unsigned char> rawColors(w*h*4);
	// for (int i = 0; i < data.segmentCount; i++){
	// 	segColorsR[i] += v[i] * 50 * (i % 12 == 0);
	// 	segColorsG[i] += v[i] * 50 * (i % 12 == 1);
	// 	segColorsB[i] += v[i] * 50 * (i % 12 == 2);
	// 	segColorsR[i] += v[i] * 100 * (i % 12 == 3);
	// 	segColorsG[i] += v[i] * 100 * (i % 12 == 4);
	// 	segColorsB[i] += v[i] * 100 * (i % 12 == 5);
	// 	segColorsR[i] += v[i] * 150 * (i % 12 == 6);
	// 	segColorsG[i] += v[i] * 150 * (i % 12 == 7);
	// 	segColorsB[i] += v[i] * 150 * (i % 12 == 8);
	// }
	printf("Number of segments: %d \n", data.segmentCount);
	for (int i = 0; i < data.segmentCount; i++){
        segColorsR[i] = rand();
        segColorsG[i] = rand();
        segColorsB[i] = rand();
    }
	for (int i = 0; i < w * h; i++) {
		//if ((i % w != 0 && data.segment[i - 1] != data.segment[i]) || (i % w + 1 != w && data.segment[i + 1] != data.segment[i])
		//	|| (i / w != 0 && data.segment[i - w] != data.segment[i]) || (i / w + 1 != h && data.segment[i + w] != data.segment[i])) {
		//	rawColors[i * 4] = 0;
		//	rawColors[i * 4 + 1] = 0;
		//	rawColors[i * 4 + 2] = 0;
		//	rawColors[i * 4 + 3] = 255;
		//}
		//else {
			rawColors[i * 4] = std::max(0, segColorsR[data.segment[i]] / v[data.segment[i]]);
			rawColors[i * 4 + 1] = std::max(0, segColorsG[data.segment[i]] / v[data.segment[i]]);
			rawColors[i * 4 + 2] = std::max(0, segColorsB[data.segment[i]] / v[data.segment[i]]);
			rawColors[i * 4 + 3] = 255;
		//}
	}
	lodepng::encode("img8_seq.png", rawColors, w, h);
}
