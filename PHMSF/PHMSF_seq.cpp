
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
	unsigned int x;
	unsigned int y;
};

struct Edge{
	unsigned int weight;
	//struct Node n1;
	//struct Node n2;
	unsigned int n1;
	unsigned int n2;
};

struct Graph{
	struct Edge* edges;
	unsigned int num_edges;
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

	int i, j;
	//can take num edges now
	//Later can move it curr_edges and then sort
	//checking only top and right edges
	int edge_cnt = 0;
	for(i = 0 ; i < width; i++){
		for(j = 0; j < height; j++){
			int curr_idx = j * width + i;

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

void swap(struct Edge* e1, struct Edge* e2){
	
	struct Edge temp;

	temp.n1 = e1->n1;
	temp.n2 = e1->n2;
	temp.weight = e1->weight;

	e1->n1 = e2->n1;
	e1->n2 = e2->n2;
	e1->weight = e2->weight;

	e2->n1 = temp.n1;
	e2->n2 = temp.n2;
	e2->weight = temp.weight;

}

void quickSort(struct Edge* edges, int l, int r)
{
	//printf("l:%d, r:%d\n", l, r);
    // Base case: No need to sort arrays of length <= 1
    if (l >= r)
    {
        return;
    }
    
    // Choose pivot to be the last element in the subarray
    int pivot = edges[r].weight;

    // Index indicating the "split" between elements smaller than pivot and 
    // elements greater than pivot
    int cnt = l;

    // Traverse through array from l to r
    for (int i = l; i <= r; i++)
    {
        // If an element less than or equal to the pivot is found...
        if (edges[i].weight <= pivot)
        {
            // Then swap arr[cnt] and arr[i] so that the smaller element arr[i] 
            // is to the left of all elements greater than pivot
            swap(&edges[cnt], &edges[i]);

            // Make sure to increment cnt so we can keep track of what to swap
            // arr[i] with
            cnt++;
        }
    }
    
    // NOTE: cnt is currently at one plus the pivot's index 
    // (Hence, the cnt-2 when recursively sorting the left side of pivot)
    quickSort(edges, l, cnt-2); // Recursively sort the left side of pivot
    quickSort(edges, cnt, r);   // Recursively sort the right side of pivot
}

bool operator<(const struct Edge &a, const struct Edge &b) {
  return a.weight < b.weight;
}

//Using quick sort
struct Graph* sort_graph(struct Graph* graph){

	//Using quick sort
	//quickSort(graph->edges, 0, (int)(graph->num_edges - 1));
	int num_edges = graph->num_edges;
	Edge *edges =  graph->edges;
	std::sort(edges, edges + num_edges);

	return graph;
}


// lower_bound_combine(struct Graph* graph, unsigned int min_weight, unsigned int max_weight){

// }

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
	printf("after sort\n");
	//combine the ones below the minimum weight
	//graph = lower_bound_combine(graph);
	//Use credit to expand the remaining regions
	//edge_heursitc();

	//Finally write to the output image
	//outImage
	

}
