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

#ifndef SEGMENT_IMAGE
#define SEGMENT_IMAGE

#include "filter.h"
#include "image.h"
#include "misc.h"
#include "segment-graph.h"
#include <cstdlib>

// random color
rgb random_rgb() {
  rgb c;
  double r;

  c.r = (uchar)rand();
  c.g = (uchar)rand();
  c.b = (uchar)rand();

  return c;
}

// dissimilarity measure between pixels
static inline float diff(image<float> *r, image<float> *g, image<float> *b,
                         int x1, int y1, int x2, int y2) {
  return sqrt(square(imRef(r, x1, y1) - imRef(r, x2, y2)) +
              square(imRef(g, x1, y1) - imRef(g, x2, y2)) +
              square(imRef(b, x1, y1) - imRef(b, x2, y2)));
}

/*
 * Segment an image
 *
 * Returns a color image representing the segmentation.
 *
 * im: image to segment.
 * sigma: to smooth the image.
 * c: constant for treshold function.
 * min_size: minimum component size (enforced by post-processing stage).
 * num_ccs: number of connected components in the segmentation.
 */
image<rgb> *segment_image(image<rgb> *im, float sigma, float c, int min_size,
                          int *num_ccs) {
  int width = im->width();
  int height = im->height();
  int pixels = width * height;

  image<float> *r = new image<float>(width, height);
  image<float> *g = new image<float>(width, height);
  image<float> *b = new image<float>(width, height);

  // smooth each color channel
  double start6 = currentSeconds();
  image<float> *smooth_r;
  image<float> *smooth_g;
  image<float> *smooth_b;
#pragma omp parallel
  {
#pragma omp for collapse(2) nowait
    for (int y = 0; y < height; y++) {
      for (int x = 0; x < width; x++) {
        imRef(r, x, y) = imRef(im, x, y).r;
        imRef(g, x, y) = imRef(im, x, y).g;
        imRef(b, x, y) = imRef(im, x, y).b;
      }
    }
    if (omp_get_thread_num() == 0)
      smooth_r = smooth(r, sigma);
    if (omp_get_thread_num() == 1)
      smooth_g = smooth(g, sigma);
    if (omp_get_thread_num() == 2)
      smooth_b = smooth(b, sigma);
  }
  double delta6 = currentSeconds() - start6;
  printf("Time taken for smoothing = %lf secs\n", delta6);
  delete r;
  delete g;
  delete b;

  // build graph
  double start5 = currentSeconds();
  edge *edges = new edge[pixels * 4];
  int num = 0;
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      if (x < width - 1) {
        edges[num].a = y * width + x;
        edges[num].b = y * width + (x + 1);
        edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x + 1, y);
        num++;
      }

      if (y < height - 1) {
        edges[num].a = y * width + x;
        edges[num].b = y * width + width + x;
        edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x, y + 1);
        num++;
      }

      if ((x < width - 1) && (y < height - 1)) {
        edges[num].a = y * width + x;
        edges[num].b = y * width + width + (x + 1);
        edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x + 1, y + 1);
        num++;
      }

      if ((x < width - 1) && (y > 0)) {
        edges[num].a = y * width + x;
        edges[num].b = y * width - width + (x + 1);
        edges[num].w = diff(smooth_r, smooth_g, smooth_b, x, y, x + 1, y - 1);
        num++;
      }
    }
  }
  double delta5 = currentSeconds() - start5;
  printf("Time taken for building graphs = %lf secs\n", delta5);
  delete smooth_r;
  delete smooth_g;
  delete smooth_b;

  // segment
  double start1 = currentSeconds();
  universe *u = segment_graph(pixels, num, edges, c);
  double delta1 = currentSeconds() - start1;
  printf("Time taken by segment_graph = %lf secs\n", delta1);

  // post process small components
  double start4 = currentSeconds();
  for (int i = 0; i < num; i++) {
    int a = u->find(edges[i].a);
    int b = u->find(edges[i].b);
    if ((a != b) && ((u->size(a) < min_size) || (u->size(b) < min_size)))
      u->join(a, b);
  }
  double delta4 = currentSeconds() - start4;
  printf("Time taken for processing small elements = %lf secs\n", delta4);
  delete[] edges;
  *num_ccs = u->num_sets();

  image<rgb> *output = new image<rgb>(width, height);

  // pick random colors for each component
  rgb *colors = new rgb[pixels];
  double start2 = currentSeconds();
  for (int i = 0; i < pixels; i++) {
    colors[i] = random_rgb();
  }
  double delta2 = currentSeconds() - start2;
  printf("Time taken by rgb_colors = %lf secs\n", delta2);

  double start3 = currentSeconds();
  for (int y = 0; y < height; y++) {
    for (int x = 0; x < width; x++) {
      int comp = u->find(y * width + x);
      imRef(output, x, y) = colors[comp];
    }
  }
  double delta3 = currentSeconds() - start3;
  printf("Time taken for copying to output = %lf secs\n", delta3);

  delete[] colors;
  delete u;

  return output;
}

#endif
