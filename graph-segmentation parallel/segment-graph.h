/*
Copyright (C) 2006 Pedro Felzenszwalb

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
*/

#ifndef SEGMENT_GRAPH
#define SEGMENT_GRAPH

#include "cycletimer.h"
#include "disjoint-set.h"
#include <algorithm>
#include <cmath>
#include <parallel/algorithm>

// threshold function
#define THRESHOLD(size, c) (c / size)

typedef struct {
  float w;
  int a, b;
} edge;

bool operator<(const edge &a, const edge &b) { return a.w < b.w; }

/*
 * Segment a graph
 *
 * Returns a disjoint-set forest representing the segmentation.
 *
 * num_vertices: number of vertices in graph.
 * num_edges: number of edges in graph
 * edges: array of edges.
 * c: constant for treshold function.
 */
universe *segment_graph(int num_vertices, int num_edges, edge *edges, float c) {
  // sort edges by weight
  double start3 = currentSeconds();
  __gnu_parallel::sort(edges, edges + num_edges);
  double delta3 = currentSeconds() - start3;
  printf("Time taken to sort edges = %lf secs\n", delta3);

  // make a disjoint-set forest
  double start1 = currentSeconds();
  universe *u = new universe(num_vertices);
  double delta1 = currentSeconds() - start1;
  printf("Time taken by universe = %lf secs\n", delta1);

  // init thresholds
  float *threshold = new float[num_vertices];
  double start = currentSeconds();
#pragma omp parallel for schedule(static, num_vertices / 8)
  for (int i = 0; i < num_vertices; i++)
    threshold[i] = THRESHOLD(1, c);
  double delta = currentSeconds() - start;
  printf("Time taken by threshold = %lf secs\n", delta);

  // for each edge, in non-decreasing weight order...
  double start2 = currentSeconds();
  for (int i = 0; i < num_edges; i++) {
    edge *pedge = &edges[i];

    // components conected by this edge
    int a = u->find(pedge->a);
    int b = u->find(pedge->b);
    if (a != b) {
      if ((pedge->w <= threshold[a]) && (pedge->w <= threshold[b])) {
        u->join(a, b);
        a = u->find(a);
        threshold[a] = pedge->w + THRESHOLD(u->size(a), c);
      }
    }
  }
  double delta2 = currentSeconds() - start2;
  printf("Time taken for segmentation algorithm loop = %lf secs\n", delta2);

  // free up
  delete threshold;
  return u;
}

#endif
