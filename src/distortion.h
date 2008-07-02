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

#ifndef __DISTORTION_H__
#define __DISTORTION_H__

#include "wcs.h"

// TODO: This is all two-dimensional.  Should be made
// multi-dimensional in the future.
#define NAXES 2

struct distortion_lookup_t {
  unsigned int  naxis[NAXES];   // size of distortion image
  double        crpix[NAXES];
  double        crval[NAXES];
  double        cdelt[NAXES];
  double       *data;
};

void
distortion_lookup_t_init(struct distortion_lookup_t* lookup);

void
distortion_lookup_t_free(struct distortion_lookup_t* lookup);

struct distortion_t {
  /* Step I:
     The CPDISja/DPja pre-linear tranformation distortion
     matrices.  May be NULL, in which case the distortion
     won't be performed.
  */
  struct distortion_lookup_t* pre_dist[NAXES];

  /* Step II:
     The CRPIXj and PCi_j or CD_i_j part of the tranformation.  If
     using CD rather than PC, pc should be set to the contents of CD,
     and CDELT should be set to 1's.
  */
  double crpix[NAXES];
  double *pc;
  int has_pc;

  /* Step III:
     The CQDISja/DQja post-linear transformation distortion matrices.
     May be NULL, in which case the distortion won't be performed.
  */
  struct distortion_lookup_t* post_dist[NAXES];

  /* Step IV:
     Recale to physical units
  */
  double cdelt[NAXES];

  /* Step V:
     CTYPE coordinate computation per agreement, handled by wcslib
  */
  struct wcsprm wcs;
};

/**
 * Initialize a distortion_t object, that describes an
 * entire transformation pipeline of Figure 1 in Paper IV.
 *
 * (Only lookup table distortions are supported, not polynomial ones.)
 *
 * TODO: docstring
 */
void
distortion_t_init(
    struct distortion_t* dist);

void
distortion_t_free(
    struct distortion_t* dist);

double
get_distortion_offset(
    const struct distortion_lookup_t *lookup,
    const double img[NAXES]);

void
distortion_pipeline(
    const struct distortion_t *dist,
    const unsigned int nelem,
    const double *pix /* [NAXES][nelem] */,
    /* Output parameters */
    double *world /* [NAXES][nelem] */);

#endif /* __DISTORTION_H__ */
