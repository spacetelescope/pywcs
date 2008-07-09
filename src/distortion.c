/*
Copyright (C) 2008 Association of Universities for Research in Astronomy (AURA)

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    1. Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    2. Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    3. The name of AURA and its representatives may not be used to
      endorse or promote products derived from this software without
      specific prior written permission.

THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

/*
 Author: Michael Droettboom
         mdroe@stsci.edu
*/

#include "distortion.h"

#include <assert.h>
#include <math.h>
#include <stdlib.h>

void
distortion_lookup_t_init(struct distortion_lookup_t* lookup) {
  unsigned int i;

  for (i = 0; i < NAXES; ++i) {
    lookup->naxis[i] = 0;
    lookup->crpix[i] = 0.0;
    lookup->crval[i] = 0.0;
    lookup->cdelt[i] = 1.0;
  }

  lookup->data = NULL;
}

void
distortion_lookup_t_free(struct distortion_lookup_t* lookup) {

}

int
distortion_t_init(struct distortion_t* dist) {
  unsigned int i;
  int status;

  assert(dist);

  for (i = 0; i < NAXES; ++i) {
    dist->pre_dist[i] = NULL;
  }

  assert(crpix);
  for (i = 0; i < NAXES; ++i) {
    dist->crpix[i] = 0.0;
  }

  dist->pc[0] = 1.0;
  dist->pc[1] = 0.0;
  dist->pc[2] = 0.0;
  dist->pc[3] = 1.0;

  for (i = 0; i < NAXES; ++i) {
    dist->post_dist[i] = NULL;
  }

  for (i = 0; i < NAXES; ++i) {
    dist->cdelt[i] = 1.0;
  }

  dist->has_pc = 0;

  status = wcsini(1, NAXES, &dist->wcs);
  assert(dist->wcs->cdelt[0] == 1.0);
  assert(dist->wcs->cdelt[1] == 1.0);
  assert(dist->wcs->crpix[0] == 0.0);
  assert(dist->wcs->crpix[1] == 0.0);
  assert(dist->wcs->pc[0]    == 1.0);
  return status;
}

int
distortion_t_free(
    struct distortion_t* dist) {
  return wcsfree(&dist->wcs);
}

/**
 * Get a value at a specific integral location in the lookup table.
 * (This is nothing more special than an array lookup with range
 * checking.)
 */
static inline double
get_dist(
    const struct distortion_lookup_t *lookup,
    const unsigned int coord[NAXES]) {
  size_t i;
  for (i = 0; i < NAXES; ++i)
    assert(coord[i] < lookup->naxis[i]);

  return *(lookup->data + (lookup->naxis[0] * coord[1]) + coord[0]);
}

/**
 * Converts a pixel coordinate to a fractional coordinate in the
 * lookup table on a single axis
 */
static inline double
image_coord_to_distortion_coord(
    const struct distortion_lookup_t *lookup,
    const unsigned int axis,
    const double img) {
  double result;

  assert(lookup);
  assert(axis < NAXES);

  result = (
      ((img - lookup->crval[axis]) / lookup->cdelt[axis]) +
      lookup->crpix[axis]);

  /* TODO: Should we just fail here? */
  if (result < 0.0)
    result = 0.0;
  else if (result >= lookup->naxis[axis])
    result = lookup->naxis[axis] - 1.0;

  return result;
}

/**
 * Converts a pixel coordinate to a fractional coordinate in the
 * lookup table.
 */
static inline void
image_coords_to_distortion_coords(
    const struct distortion_lookup_t *lookup,
    const double img[NAXES],
    /* Output parameters */
    double dist[NAXES]) {
  size_t i;

  assert(lookup);
  assert(img);
  assert(dist);

  for (i = 0; i < NAXES; ++i) {
    dist[i] = image_coord_to_distortion_coord(lookup, i, img[i]);
  }
}

/**
 * Helper function for get_distortion_offset
 */
static inline double
calculate_weight(double iw, unsigned int i0, double jw, unsigned int j0) {
  if (!i0)
    iw = 1.0 - iw;
  if (!j0)
    jw = 1.0 - jw;
  return iw * jw;
}

double
get_distortion_offset(
    const struct distortion_lookup_t *lookup,
    const double img[NAXES]) {
  double       dist[NAXES];
  double       dist_floor[NAXES];
  unsigned int dist_i[NAXES];
  unsigned int coord[NAXES];
  double       dist_weight[NAXES];
  double       result;
  size_t       i, k, l;

  assert(lookup);
  assert(img);

  image_coords_to_distortion_coords(lookup, img, dist);

  for (i = 0; i < NAXES; ++i) {
    dist_floor[i] = floor(dist[i]);
    dist_i[i] = (unsigned int)dist_floor[i];
    dist_weight[i] = dist[i] - dist_floor[i];
  }

  result = 0.0;
  for (k = 0; k < 2; ++k) {
    for (l = 0; l < 2; ++l) {
      coord[0] = dist_i[0] + l;
      coord[1] = dist_i[1] + k;
      result += (get_dist(lookup, coord) *
                 calculate_weight(dist_weight[0], l, dist_weight[1], k));
    }
  }

  return result;
}

int
distortion_pipeline(
    const struct distortion_t *dist,
    const unsigned int nelem,
    const double *pix /* [NAXES][nelem] */,
    /* Output parameters */
    double *world /* [NAXES][nelem] */) {
  const double *xi;
  const double *yi;
  double       *xo, *yo;
  double        offset;
  double        pixcrd[NAXES];
  double        pre_dist;
  double        inter[NAXES];
  double        post_dist;
  double        scaled[NAXES];
  double        worldcrd[NAXES];
  int           stat[NAXES];
  int           status;
  size_t        i, j, k;

  assert(dist);
  assert(pix);
  assert(world);

  if (dist->pc[0] == 0.0 || dist->pc[3] == 0.0) {
    return 3; // linear transformation matrix is singular
  }

  if (dist->cdelt[0] == 0.0 || dist->cdelt[1] == 0.0) {
    return 6; // Invalid coordinate transformation parameters
  }

  xi = pix;
  yi = pix + nelem;

  xo = world;
  yo = world + nelem;

  for (k = 0; k < nelem; ++k, ++xi, ++yi, ++xo, ++yo) {
    pixcrd[0] = *xi;
    pixcrd[1] = *yi;

    for (i = 0; i < NAXES; ++i) {
      /* STEP I: CPDIS distortion correction */
      if (dist->pre_dist[i]) {
        offset = get_distortion_offset(dist->pre_dist[i], pixcrd);
        pre_dist = pixcrd[i] + offset;
      } else {
        pre_dist = pixcrd[i];
      }

      /* STEP II: CRPIX linear transformation */
      inter[i] = 0.0;
      offset = pre_dist - dist->crpix[i];
      for (j = 0; j < NAXES; ++j) {
        inter[i] += dist->pc[(i<<1)+j] * offset;
      }
    }

    /* STEP III: CQDIS distortion correction */
    for (i = 0; i < NAXES; ++i) {
      if (dist->post_dist[i]) {
        offset = get_distortion_offset(dist->post_dist[i], inter);
        post_dist = inter[i] + offset;
      } else {
        post_dist = inter[i];
      }

      if (dist->has_pc) {
        /* STEP IV: CDELT scaling */
        scaled[i] = post_dist * dist->cdelt[i];
      }
    }

    /* STEP V: CTYPE coordinate computation per agreement */
    /* This is the most complicated part -- so we defer to wcslib */
    /* TODO: Error handling */
    status = wcsp2s((struct wcsprm *)&dist->wcs, NAXES, 1,
                    scaled, pixcrd, pixcrd, pixcrd, worldcrd, stat);
    if (status != 0 && status != 8)
      return status;

    *xo = worldcrd[0];
    *yo = worldcrd[1];
  }

  return 0;
}
