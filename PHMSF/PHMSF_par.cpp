
#include <string.h>
#include <getopt.h>

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <stdarg.h>
#include <ctype.h>
#include <math.h>
#include <algorithm>
#include <parallel/algorithm>

#include "lodepng.h"
#include "cycletimer.h"

#define outImg_extension "_par"
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
	//unsigned int linear_idx;
	unsigned int rank;
	unsigned int parent;
	//unsigned int x;
	//unsigned int y;
};

struct Edge{
	int weight;
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

	double start = currentSeconds();
	lodepng::decode(tmp, w, h, inpImg);
	double delta = currentSeconds() - start;
	printf("Time taken for loadepng:: decode = %lf secs\n", delta);
	struct Image img;
	img.width = w;
	img.height = h;

	img.pixels = new Pixel[w*h];
	#pragma omp parallel for schedule (static)
	for (int i = 0; i < w*h; i++){
		img.pixels[i].R = tmp[i * 4];
		img.pixels[i].G = tmp[i * 4 + 1];
		img.pixels[i].B = tmp[i * 4 + 2];
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

int getL2norm(struct Pixel p1, struct Pixel p2){
	int rDiff = (int)(p1.R - p2.R);
	int gDiff = (int)(p1.G - p2.G);
	int bDiff = (int)(p1.B - p2.B);
	rDiff = rDiff*rDiff;
	gDiff = gDiff*gDiff;
	bDiff = bDiff*bDiff;
	int res = 4 *sqrt(rDiff + gDiff + bDiff);
	return res;
}

struct Graph* build_graph(struct Image image, int maxWeight){
	struct Graph* graph = new Graph(); 
	unsigned int height = image.height;
	unsigned int width = image.width;
	int num_edges = width * height * 2 - width - height;
	graph->edges = new Edge[num_edges];
	graph->num_nodes = width* height;
	graph->nodes = new Node[graph->num_nodes];


	//use tiling to build graph
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		int num_threads = omp_get_num_threads();
		int edge_cnt;
		//Splitting into horizontal graph
		edge_cnt = tid * (height*(2*width -1)/num_threads) - width;
		int istart = tid*(height/num_threads);
		if(tid == 0){
			istart = 0;
			edge_cnt = 0;
		}
		int iend = (tid+1)*(height/num_threads);
 
	//#pragma omp parallel for
	for(int i = istart ; i < iend; i++){
		for(int j = 0; j < width; j++){
			int curr_idx = i * width + j;

			graph->nodes[curr_idx].rank = 1;
			graph->nodes[curr_idx].parent = curr_idx;

			//Check the edge case
			if(i + 1 < height){
				int top_idx = curr_idx + width;

				int weight = getL2norm(image.pixels[curr_idx], image.pixels[top_idx]);
				if(weight <= maxWeight){
					graph->edges[edge_cnt].n1 = curr_idx;
					graph->edges[edge_cnt].n2 = top_idx; 
					graph->edges[edge_cnt].weight = weight;
					edge_cnt++;
				}
			}
			if(j + 1 < width){
				int right_idx = curr_idx + 1;

				int weight = getL2norm(image.pixels[curr_idx], image.pixels[right_idx]);
				if(weight <= maxWeight){
					graph->edges[edge_cnt].n1 = curr_idx;
					graph->edges[edge_cnt].n2 = right_idx; 
					graph->edges[edge_cnt].weight = weight;
					edge_cnt++;
				}
			}
		}
	}
	int maxEdges = (tid+1) * (height*(2*width -1)/num_threads) - width;
	//printf("val at %d\n", graph->edges[edge_cnt-1].weight);
	//if(edge_cnt < maxEdges){
		//printf("Differnce %d\n", maxEdges - edge_cnt);
		//memset(graph->edges + num_edges, 0, sizeof(struct Edge)*(maxEdges-edge_cnt));
	//}
	//printf("after val at %d\n", graph->edges[edge_cnt-1].weight);
	if(tid == num_threads - 1){
		printf("NUMBER OF EDGES = %d\n", edge_cnt);
		num_edges = edge_cnt;
	}
	}
	//printf("NUMBER OF EDGES = %d\n", edge_cnt);
	graph->num_edges = num_edges;
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
	return n1;
}
#define computeCredit(regionSize) sqrt(4 * 3.14 * regionSize)

struct Region* find_regions(unsigned int width, 
	unsigned int height, unsigned int minRegionSize, struct Node* nodes){
	int size = width * height;

	struct Region* regions = new Region[width*height];
	printf("number of threads: %d\n", omp_get_num_threads());

	#pragma omp parallel for schedule (static)
	for (int i = 0; i < width*height; i++) {
		regions[i].size = regions[i].rep = regions[i].credit = 0;
	}
	
	#pragma omp parallel for schedule (static)
	for (int j = 0; j < width*height; j++) {
		int parent = union_find(nodes, j);
		regions[parent].size++;
		regions[j].rep = parent;
	}

	#pragma omp parallel for schedule (static)
	for (int i = 0; i < width*height; i++) {
		if (regions[i].size >= minRegionSize)
			regions[i].credit = computeCredit(regions[i].size);
		else{

			regions[i].credit = sqrt(width*height)*20; 
		}
	}
	return regions;

}

bool operator<(const struct Edge &a, const struct Edge &b) {
  return a.weight < b.weight;
}

//Using quick sort
struct Graph* sort_graph(struct Graph* graph){

	int num_edges = graph->num_edges;
	Edge *edges =  graph->edges;
	__gnu_parallel::sort(edges, edges + num_edges);
	//std::sort(edges, edges + num_edges);
	return graph;
}


struct Graph* lower_bound_combine(struct Graph* graph, 
	unsigned int min_weight){
	int i;
	int united_cnt = 0;
	#pragma omp parallel for schedule (static)
	for(i = 0; i < graph->num_edges; i++){
		if(graph->edges[i].weight <= min_weight){

			union_unite(graph->nodes, graph->edges[i].n1, graph->edges[i].n2);
		}
	}
	return graph;
 }



void edge_heursitc(struct Graph* graph, struct Region* regions, unsigned int minWeight){

 	struct Edge* edges = graph->edges;
 	unsigned int edge_cnt = graph->num_edges;
 	printf("edge cnt: %d\n", edge_cnt);
 	//int i;
 	//#pragma omp parallel for schedule (static)
 	//Need to figure out all the border edges
 #pragma omp parallel
 { 
 		int tid = omp_get_thread_num();
		int num_threads = omp_get_num_threads();
		//Splitting into horizontal graph
		int istart = tid*(edge_cnt/num_threads); 
		int iend = std::min((tid+1)*(edge_cnt/num_threads), edge_cnt); 
	#pragma omp parallel for
 	for(int i = istart; i < iend; i++){
 		if(edges[i].weight <= minWeight){
 			continue;
 		}
 		unsigned int a = union_find(graph->nodes, edges[i].n1);
		unsigned int b = union_find(graph->nodes, edges[i].n2);
 		if(a != b){
			int credit = std::min(regions[a].credit, regions[b].credit);
			if (credit > edges[i].weight) {
				int s = union_unite(graph->nodes, a, b);
				#pragma omp atomic
				regions[s].credit = credit - edges[i].weight;
			}

 		}
 	}
 }
 }

void smooth_image(struct Pixel* pixel, unsigned int width, unsigned int height){
 	int i, j, k, l;
	for(i = 2 ; i < height-2; i++){
		for(j = 2; j < width-2; j++){
			int curr_idx = i * width + j;

			int tempR = 0;//pixel[curr_idx].R;
			int tempG = 0;//pixel[curr_idx].G;
			int tempB = 0;//pixel[curr_idx].B;
			for(k = -2; k < 3 ; k++){
					tempR += pixel[curr_idx + k].R;
					tempG += pixel[curr_idx + k].G;
					tempB += pixel[curr_idx + k].B;  
			}
			for(k = -2; k < 3 ; k++){
				//if(k == 0) continue;
					tempR += pixel[curr_idx + k*width].R;
					tempG += pixel[curr_idx + k*width].G;
					tempB += pixel[curr_idx + k*width].B;  
			}
			pixel[curr_idx].R = tempR/10;
			pixel[curr_idx].G = tempG/10;
			pixel[curr_idx].B = tempB/10;
		}

	}
}

void output_image(struct Graph *graph, unsigned int width, unsigned int height, std::string outImg){

	struct Pixel* colours = new Pixel[width*height];
	int num_segments = 0;
#pragma omp parallel
{
	unsigned int threadSeed = omp_get_thread_num() * 10;

	#pragma omp for schedule(static)
	for (int i = 0; i < width*height; i++) {
		if (i == union_find(graph->nodes, i)) {
			colours[i].R = rand_r(&threadSeed);
			colours[i].G = rand_r(&threadSeed);
			colours[i].B = rand_r(&threadSeed);
			#pragma omp atomic
			num_segments++;
		}
	}
}
	printf("Number of Segments: %d\n", num_segments);

	std::vector<unsigned char> rawColors(width*height*4);

	#pragma omp for schedule(static)
	for (int i = 0; i < width*height; i++) {

		int parent = union_find(graph->nodes, i);

		rawColors[i * 4] = colours[parent].R;
		rawColors[i * 4 + 1] = colours[parent].G;
		rawColors[i * 4 + 2] = colours[parent].B;
		rawColors[i * 4 + 3] = 255;
	}
	lodepng::encode(outImg, rawColors, width, height);
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

	double start, delta;

	char* optstring = "hg:i:m:M:r:s"; 
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
				start = currentSeconds();
				image = loadImage(inpImg);
				delta = currentSeconds()-start;
				printf("Time taken for loading image = %lf secs\n", delta);

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
			case 's':
				omp_set_num_threads(1);
		}
	}
	if(inpFileExists == false){
		printf("Please provide an input image\n");
		usage();
	}

	
	/*add gaussian noise*/
	start = currentSeconds();
	smooth_image(image.pixels, image.width, image.height);
	delta = currentSeconds()-start;
	printf("Time taken for smoothening image = %lf secs\n", delta);
	
	//add graph
	start = currentSeconds();
	graph = build_graph(image, maxWeight);
	delta = currentSeconds()-start;
	printf("Time taken for building graph = %lf secs\n", delta);


	//sort edges according to their weights
	start = currentSeconds();
	graph = sort_graph(graph);
	delta = currentSeconds()-start;
	printf("Time taken for sorting graph = %lf secs\n", delta);

	start = currentSeconds();
	graph = lower_bound_combine(graph, minWeight);
	delta = currentSeconds()-start;
	printf("Time taken for combining regions using lower bounds = %lf secs\n", delta);
	
	//combine the ones below the minimum weight
	start = currentSeconds();
	struct Region* regions = find_regions(image.width, image.height, minRegion, graph->nodes);
	delta = currentSeconds()-start;
	printf("Time taken for find regions = %lf secs\n", delta);
	
	//Use credit to expand the remaining regions
	start = currentSeconds();
	edge_heursitc(graph, regions, minWeight);
	delta = currentSeconds()-start;
	printf("Time taken for edge heursitc = %lf secs\n", delta);

	//Write to the output image
	start = currentSeconds();
	output_image(graph, image.width, image.height, outImg);
	delta = currentSeconds()-start;
	printf("Time taken for img output = %lf secs\n", delta);
}
