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

/* TODO: n-dimensional support */

int
distortion_lookup_t_init( distortion_lookup_t* lookup) {
  unsigned int i;

  for (i = 0; i < NAXES; ++i) {
    lookup->naxis[i] = 0;
    lookup->crpix[i] = 0.0;
    lookup->crval[i] = 0.0;
    lookup->cdelt[i] = 1.0;
  }

  lookup->data = NULL;

  return 0;
}

int
distortion_lookup_t_free( distortion_lookup_t* lookup) {
  return 0;
}

/**
 * Get a value at a specific integral location in the lookup table.
 * (This is nothing more special than an array lookup with range
 * checking.)
 */
static inline float
get_dist(
    const  distortion_lookup_t *lookup,
    const unsigned int coord[NAXES]) {
  unsigned int cropped[NAXES];
  unsigned int i;

  for (i = 0; i < NAXES; ++i) {
    cropped[i] = coord[i] >= lookup->naxis[i] ? lookup->naxis[i] - 1 : coord[i];
  }

  return *(lookup->data + (lookup->naxis[0] * cropped[1]) + cropped[0]);
}

/**
 * Converts a pixel coordinate to a fractional coordinate in the
 * lookup table on a single axis
 */
static inline double
image_coord_to_distortion_coord(
    const  distortion_lookup_t *lookup,
    const unsigned int axis,
    const double img) {
  double result;

  assert(lookup);
  assert(axis < NAXES);

  result = (
      ((img - lookup->crval[axis]) / lookup->cdelt[axis]) +
      lookup->crpix[axis]);

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
    const  distortion_lookup_t *lookup,
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
  assert(iw >= 0.0 && iw < 1.0);
  assert(jw >= 0.0 && jw < 1.0);

  if (!i0) iw = 1.0 - iw;
  if (!j0) jw = 1.0 - jw;
  return iw * jw;
}

double
get_distortion_offset(
    const distortion_lookup_t *lookup,
    const double img[NAXES]) {
  double       dist[NAXES];
  double       dist_floor[NAXES];
  unsigned int dist_ifloor[NAXES];
  unsigned int coord[NAXES];
  double       dist_weight[NAXES];
  double       result;
  size_t       i, k, l;

  assert(lookup);
  assert(img);

  image_coords_to_distortion_coords(lookup, img, dist);

  for (i = 0; i < NAXES; ++i) {
    dist_floor[i] = floor(dist[i]);
    dist_ifloor[i] = (unsigned int)dist_floor[i];
    dist_weight[i] = dist[i] - dist_floor[i];
  }

  result = 0.0;
  for (k = 0; k < 2; ++k) {
    for (l = 0; l < 2; ++l) {
      coord[0] = dist_ifloor[0] + l;
      coord[1] = dist_ifloor[1] + k;
      result += (double)get_dist(lookup, coord) * \
                calculate_weight(dist_weight[0], l, dist_weight[1], k);
    }
  }

  return result;
}

int
p4_pix2foc(
    const unsigned int naxes,
    const distortion_lookup_t *lookup[NAXES], /* [NAXES] */
    const unsigned int nelem,
    const double* pix, /* [NAXES][nelem] */
    double *foc /* [NAXES][nelem] */) {
  unsigned int i, j;

  assert(naxes == NAXES);
  assert(lookup);
#ifndef NDEBUG
  for (i = 0; i < naxes; ++i) {
    assert(lookup[i]);
    assert(lookup[i]->data);
  }
#endif
  assert(pix);
  assert(foc);
  assert(pix != foc);

  if (pix == NULL || foc == NULL || lookup[0] == NULL || lookup[1] == NULL) {
    return 1;
  }

  for (j = 0; j < nelem; ++j) {
    for (i = 0; i < NAXES; ++i) {
      foc[i] = pix[i] + get_distortion_offset(lookup[i], pix);
    }
    pix += naxes;
    foc += naxes;
  }

  return 0;
}

