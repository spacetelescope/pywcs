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
#include <string.h>

/* TODO: n-dimensional support */

int
distortion_lookup_t_init(
    distortion_lookup_t* lookup) {

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

void
distortion_lookup_t_free(
    /*@unused@*/ distortion_lookup_t* lookup) {

  /*@empty@*/
}

/**
 * Get a value at a specific integral location in the lookup table.
 * (This is nothing more special than an array lookup with range
 * checking.)
 */
static inline float
get_dist_clamp(
    const distortion_lookup_t * const lookup,
    const int x,
    const int y) {

  return *(lookup->data +
           (lookup->naxis[0] * CLAMP(y, 0, lookup->naxis[1] - 1)) +
           CLAMP(x, 0, lookup->naxis[0] - 1));
}

static inline float
get_dist(
    const distortion_lookup_t * const lookup,
    const int x,
    const int y) {

  return *(lookup->data + (lookup->naxis[0] * y) + x);
}

/**
 * Converts a pixel coordinate to a fractional coordinate in the
 * lookup table on a single axis
 */
static inline double
image_coord_to_distortion_coord(
    const distortion_lookup_t * const lookup,
    const unsigned int axis,
    const double img) {

  double result;

  assert(lookup != NULL);
  assert(axis < NAXES);

  result = (
      ((img - lookup->crval[axis]) / lookup->cdelt[axis]) +
      lookup->crpix[axis]);

  return CLAMP(result, 0.0, (double)(lookup->naxis[axis] - 1));
}

/**
 * Converts a pixel coordinate to a fractional coordinate in the
 * lookup table.
 */
static inline void
image_coords_to_distortion_coords(
    const distortion_lookup_t * const lookup,
    const double * const img /* [NAXES] */,
    /* Output parameters */
    /*@out@*/ double *dist /* [NAXES] */) {

  unsigned int i;

  assert(lookup != NULL);
  assert(img != NULL);
  assert(dist != NULL);

  for (i = 0; i < NAXES; ++i) {
    dist[i] = image_coord_to_distortion_coord(lookup, i, img[i]);
  }
}

/* /\** */
/*  * Helper function for get_distortion_offset */
/*  *\/ */
/* static inline double */
/* calculate_weight( */
/*     double iw, */
/*     const unsigned int i0, */
/*     double jw, */
/*     const unsigned int j0) { */

/*   assert(iw >= 0.0 && iw < 1.0); */
/*   assert(jw >= 0.0 && jw < 1.0); */

/*   if (i0 == 0) iw = 1.0 - iw; */
/*   if (j0 == 0) jw = 1.0 - jw; */
/*   return iw * jw; */
/* } */

/* #define CALCULATE_WEIGHT(iw, i0, jw, j0) (((i0 == 0) ? (1.0 - iw) : iw) * ((j0 == 0) ? (1.0 - jw) : jw)) */

double
get_distortion_offset(
    const distortion_lookup_t * const lookup,
    const double * const img /*[NAXES]*/) {

  double       dist[NAXES];
  double       dist_floor[NAXES];
  int          dist_ifloor[NAXES];
  double       dist_weight[NAXES];
  double       dist_iweight[NAXES];
  double       result;
  unsigned int i;

  assert(lookup != NULL);
  assert(img != NULL);

  image_coords_to_distortion_coords(lookup, img, dist);

  for (i = 0; i < NAXES; ++i) {
    dist_floor[i] = floor(dist[i]);
    dist_ifloor[i] = (int)dist_floor[i];
    dist_weight[i] = dist[i] - dist_floor[i];
    dist_iweight[i] = 1.0 - dist_weight[i];
  }

  /* If we may need to clamp the lookups, use this slower approach */
  if (dist_ifloor[0] < 0 ||
      dist_ifloor[1] < 0 ||
      dist_ifloor[0] >= lookup->naxis[0] - 1 ||
      dist_ifloor[1] >= lookup->naxis[1] - 1) {
    result =
      (double)get_dist_clamp(lookup, dist_ifloor[0], dist_ifloor[1]) * dist_weight[0] * dist_weight[1] +
      (double)get_dist_clamp(lookup, dist_ifloor[0], dist_ifloor[1] + 1) * dist_weight[0] * dist_iweight[1] +
      (double)get_dist_clamp(lookup, dist_ifloor[0] + 1, dist_ifloor[1]) * dist_iweight[0] * dist_weight[1] +
      (double)get_dist_clamp(lookup, dist_ifloor[0] + 1, dist_ifloor[1] + 1) * dist_iweight[0] * dist_iweight[1];
  /* Else, we don't need to clamp 4 times for each pixel */
  } else {
    result =
      (double)get_dist(lookup, dist_ifloor[0], dist_ifloor[1]) * dist_weight[0] * dist_weight[1] +
      (double)get_dist(lookup, dist_ifloor[0], dist_ifloor[1] + 1) * dist_weight[0] * dist_iweight[1] +
      (double)get_dist(lookup, dist_ifloor[0] + 1, dist_ifloor[1]) * dist_iweight[0] * dist_weight[1] +
      (double)get_dist(lookup, dist_ifloor[0] + 1, dist_ifloor[1] + 1) * dist_iweight[0] * dist_iweight[1];
  }

  return result;
}

int
p4_pix2deltas(
    const unsigned int naxes,
    const distortion_lookup_t **lookup, /* [NAXES] */
    const unsigned int nelem,
    const double* pix, /* [NAXES][nelem] */
    double *foc /* [NAXES][nelem] */) {

  int j;

  assert(naxes == NAXES);
  assert(lookup != NULL);
  assert(pix != NULL);
  assert(foc != NULL);
#ifndef NDEBUG
  unsigned int k;
  for (k = 0; k < naxes; ++k) {
    if (lookup[k] != NULL) {
      assert(lookup[k]->data != NULL);
    }
  }
#endif

  if (pix == NULL || foc == NULL) {
    return 1;
  }

#pragma omp parallel
  {
    int i;
    double* foc0;
    const double* pix0;

#pragma omp for
    for (j = 0; j < nelem; ++j) {
      pix0 = pix + (j * NAXES);
      foc0 = foc + (j * NAXES);
      for (i = 0; i < NAXES; ++i) {
        if (lookup[i]) {
          foc0[i] += get_distortion_offset(lookup[i], pix0);
        }
      }
    }
  } /* end of parallel section  */

  return 0;
}

int
p4_pix2foc(
    const unsigned int naxes,
    const distortion_lookup_t **lookup, /* [NAXES] */
    const unsigned int nelem,
    const double* pix, /* [NAXES][nelem] */
    double *foc /* [NAXES][nelem] */) {

  assert(pix);
  assert(foc);

  if (pix != foc) {
      memcpy(foc, pix, sizeof(double) * naxes * nelem);
  }

  return p4_pix2deltas(naxes, lookup, nelem, pix, foc);
}

