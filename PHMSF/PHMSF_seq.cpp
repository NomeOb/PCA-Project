
#include <string.h>
#include <getopt.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include <algorithm>

#include "lodepng.h"

#define outImg_extension "_seq"
#define file_extension ".png"
#define DEFAULT_MIN_WEIGHT 20
#define DEFAULT_MAX_WEIGHT 90
#define DEFAULT_MIN_REGION_SIZE 50

#define L2_POWER_FACTOR 2

struct Pixel{
	unsigned char R;
	unsigned char G;
	unsigned char B;
};

struct segmentedImage{
	unsigned int num_segments;
};

struct Image{

	unsigned int width;
	unsigned int height;
	struct Pixel* pixels;
};

struct Node{
	unsigned int linear_idx;
	unsigned int rank;
	unsigned int parent;
	//unsigned int x;
	//unsigned int y;
};

struct Edge{
	unsigned int weight;
	//struct Node n1;
	//struct Node n2;
	unsigned int n1;
	unsigned int n2;
};

struct Graph{
	struct Node *nodes;
	struct Edge* edges;
	unsigned int num_edges;
	unsigned int num_nodes;
};

struct Region{
	unsigned int size;
	unsigned int rep;
	unsigned int credit;
};

void usage(){
	printf("Usage: [-h] [-i IFILE] [-m MIN] [-M MAX] [-r MINSIZE]\n");
    printf("   -h         Print this message\n");
    printf("   -i IFILE   Input Image file (.png format)\n");
    printf("   -m MIN     Minimum Edge Weight Difference\n");
    printf("   -M MAX     Maximum Edge Weight Difference\n");
    printf("   -r MINSIZE Minimum Region Size\n");
    exit(0);

}

bool ends_with(std::string a, std::string b) {
	int len = b.length();
	int pos = a.length() - len;

	std::string end = a.substr (pos,len);

	return !(b.compare(end));
}

bool valid_file(char *optarg){
	FILE* inpImg = NULL;
	std::string fileType = file_extension;

	if(ends_with(optarg, fileType)){
		inpImg = fopen(optarg, "r");
	}
	else {
		printf("Invalid File format %s, please input a .png file\n", optarg);
		return false;
	}
	if(inpImg == NULL){
		printf("Unable to find file %s\n", optarg);
		return false;
	}
	else {
		fclose(inpImg);
		return true;
	}
}

struct Image loadImage(std::string inpImg){

	std::vector<unsigned char> tmp;
	unsigned int w, h;
	lodepng::decode(tmp, w, h, inpImg);
	struct Image img;
	img.width = w;
	img.height = h;

	img.pixels = new Pixel[w*h];
	for (int i = 0; i < w*h; i++){
		img.pixels->R = tmp[i * 4];
		img.pixels->G = tmp[i * 4 + 1];
		img.pixels->B = tmp[i * 4 + 2];
	}
	return img;
}

std::string outImageName(std::string outImg){
	//.png is 4 letters
	std::string ext = file_extension;
	int insertAt = outImg.length() - ext.length();
	outImg.insert(insertAt, outImg_extension);
	return outImg;
}

unsigned int getL2norm(struct Pixel* n1, struct Pixel* n2){
	return sqrt(pow((n1->R - n2->R), L2_POWER_FACTOR) + pow((n1->G - n2->G), L2_POWER_FACTOR) + pow((n1->B - n2->B), L2_POWER_FACTOR));
}

struct Graph* build_graph(struct Image image){
	struct Graph* graph = new Graph(); //can do malloc also
	unsigned int height = image.height;
	unsigned int width = image.width;
	graph->num_edges = width * height * 2 - image.width - image.height;
	graph->edges = new Edge[graph->num_edges];
	graph->num_nodes = width* height;
	graph->nodes = new Node[graph->num_nodes];

	int i, j;
	//can take num edges now
	//Later can move it curr_edges and then sort
	//checking only top and right edges
	int edge_cnt = 0;
	for(i = 0 ; i < width; i++){
		for(j = 0; j < height; j++){
			int curr_idx = j * width + i;

			//Might not need linear index
			graph->nodes[curr_idx].linear_idx = curr_idx;
			graph->nodes[curr_idx].rank = 1;
			graph->nodes[curr_idx].parent = curr_idx;


			//Check the edge case
			if(j + 1 < height){
				int top_idx = curr_idx + width;
				graph->edges[edge_cnt].n1 = curr_idx;
				graph->edges[edge_cnt].n2 = top_idx; 
				graph->edges[edge_cnt].weight = getL2norm(&image.pixels[curr_idx], &image.pixels[top_idx]);
				edge_cnt++;
			}
			if(i + 1 < width){
				int right_idx = curr_idx + 1;
				graph->edges[edge_cnt].n1 = curr_idx;
				graph->edges[edge_cnt].n2 = right_idx; 
				graph->edges[edge_cnt].weight = getL2norm(&image.pixels[curr_idx], &image.pixels[right_idx]);
				edge_cnt++;
			}
		}
	}
	graph->num_edges = edge_cnt;
	return graph;
}


unsigned int union_find(struct Node* nodes, unsigned int node_id){
	return nodes[node_id].parent = (node_id == nodes[node_id].parent ? node_id : union_find(nodes, nodes[node_id].parent));
}

unsigned int union_unite(struct Node *nodes, unsigned int n1, unsigned int n2){
	n1 = union_find(nodes, n1);
	n2 = union_find(nodes, n2);
	if(n1 != n2){
		if(nodes[n2].rank > nodes[n1].rank){
			std::swap(n1, n2);
		}
		nodes[n2].parent = nodes[n1].parent;
		nodes[n1].rank += nodes[n2].rank;
	}
	//maybe check for same ranks here??
	return n1;
}
int computeCredit(int regionSize) {
	return sqrt(4 * 3.14 * regionSize);
}

struct Region* find_regions(unsigned int width, 
	unsigned int height, unsigned int minRegionSize, struct Node* nodes){
	int size = width * height;
	printf("Size :%d\n", size);
	struct Region* regions = new Region[width*height];

	//memset(regions, 0, sizeof(regions)*height*width);
	int i;
	for(i = 0; i < width* height ; i ++){
		regions[i].size = 0;
		regions[i].credit = 0;
		regions[i].rep = 0;
	}
	printf("Error?\n");
	for (i = 0; i < width*height; i++) {
		printf("from 1st for i %d\n", i);
		regions[union_find(nodes, i)].size++;
		regions[union_find(nodes, i)].rep = union_find(nodes, i);
	}
	printf("Error first for?\n");
	for (int i = 0; i < width*height; i++) {
		if (regions[i].size >= minRegionSize)
			regions[i].credit = computeCredit(regions[i].size);
		else
			regions[i].credit = 100000;
	}
	printf("Error first second?\n");
	return regions;

}

bool operator<(const struct Edge &a, const struct Edge &b) {
  return a.weight < b.weight;
}

//Using quick sort
struct Graph* sort_graph(struct Graph* graph){

	//Using quick sort
	//quickSort(graph->edges, 0, (int)(graph->num_edges - 1));
	//quick sort was taking too long so switched to std::sort
	int num_edges = graph->num_edges;
	Edge *edges =  graph->edges;
	std::sort(edges, edges + num_edges);

	return graph;
}


struct Graph* lower_bound_combine(struct Graph* graph, 
	unsigned int min_weight){
	int i;
	for(i = 0; i < graph->num_edges; i++){
		if(graph->edges[i].weight <= min_weight){
			union_unite(graph->nodes, graph->edges[i].n1, graph->edges[i].n2);

		}
	}

 }

bool isUnited(struct Node* nodes, unsigned int a, unsigned int b) {
		return union_find(nodes, a) == union_find(nodes, b);
	}


void edge_heursitc(struct Graph* graph, struct Region* regions, unsigned int minWeight){

 	struct Edge* edges = graph->edges;
 	unsigned int edge_cnt = graph->num_edges;
 	int i;
 	for(i = 0; i < edge_cnt; i++){
 		if(edges[i].weight >= minWeight){
 			continue;
 		}
 		if(!isUnited(graph->nodes, edges[i].n1, edges[i].n2)){
 			unsigned int a = union_find(graph->nodes, edges[i].n1);
			unsigned int b = union_find(graph->nodes, edges[i].n2);
			int credit = std::min(regions[a].credit, regions[b].credit);
			if (credit > edges[i].weight) {
				int s = union_unite(graph->nodes, a, b);
				regions[s].credit = credit - edges[i].weight;
			}

 		}
 	}

 }

int main(int argc, char *argv[]){

	bool inpFileExists = false;
	struct Image image;
	struct segmentedImage outImage;
	struct Graph *graph;

	unsigned int minWeight = DEFAULT_MIN_WEIGHT;
	unsigned int maxWeight = DEFAULT_MAX_WEIGHT;
	unsigned int minRegion = DEFAULT_MIN_REGION_SIZE;

	std::string inpImg;
	std::string outImg;

	char* optstring = "hg:i:m:M:r"; //"hg:r:R:n:s:i:q";
	char c;

	while((c = getopt(argc, argv, optstring)) != -1){
		switch(c){
			case 'h':
				usage();
				break;
			case 'i':
				//check if the file exists:
				inpFileExists = valid_file(optarg);
				if(!inpFileExists){
					exit(0);
				}
				inpImg = optarg;
				//load the PNG image
				image = loadImage(inpImg);
				outImg = outImageName(inpImg);
				printf("output:%s\n", outImg.c_str());
				printf("input:%s\n", inpImg.c_str());
				printf("width: %d, height: %d \n", image.width, image.height);
				 
			break;
			case 'm':
				minWeight = atoi(optarg);
			break;
			case 'M':
				maxWeight = atoi(optarg);
			break;
			case 'r':
				minRegion = atoi(optarg);
			break;
		}
	}
	if(inpFileExists == false){
		printf("Please provide an input image\n");
		usage();
	}
	/*add gaussian noise*/
	//add_gaussian_noise(image.pixels);
	
	//add graph
	graph = build_graph(image);
	//sort edges according to their weights
	graph = sort_graph(graph);
	//graph now has the sorted edges. 
	//add nodes to it

	//graph = mege_nodes(graph);

	graph = lower_bound_combine(graph, minWeight);
	//combine the ones below the minimum weight
	//error in this function
	printf("Lower Bound?? \n");
	//printf()
	struct Region* regions = find_regions(image.width, image.height, minRegion, graph->nodes);
	printf("find_regions?\n");	
	//Use credit to expand the remaining regions
	edge_heursitc(graph, regions, minWeight);
	printf("edges?\n");
	//count the number of regions
	//make another array??
	//Finally write to the output image
	//outImage(regions)
}
