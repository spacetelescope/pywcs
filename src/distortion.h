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

/*
 MGDTODO: It's preferable to do the pipelining of SIP + Paper IV
 distortion + WCS in Python.  In which case, it's much easier to drop
 the Distortion class (which does the whole pipeline of Paper IV,
 Figure 1), in favor of just doing the first box in that figure.

 The complication is that it is hard to insert steps into the wcslib
 processing, so CQDIS becomes harder (or impossible) to handle.  For
 now, we are assuming we don't need CQDIS.  If that really is the
 case, the whole Distortion class (but not DistortionLookupTable)
 could be removed.

 Define DISTORTION_PIPELINE to compile the old Distortion class.
*/

#ifndef __DISTORTION_H__
#define __DISTORTION_H__

#ifdef DISTORTION_PIPELINE
#include "wcs.h"
#endif

/* TODO: This is all two-dimensional.  Should be made
   multi-dimensional in the future. */
#define NAXES 2

#define MAXAXES 6

/**
 * A structure to contain the information for a single distortion lookup table
 */
struct distortion_lookup_t {
  unsigned int  naxis[NAXES];   /* size of distortion image */
  double        crpix[NAXES];
  double        crval[NAXES];
  double        cdelt[NAXES];
  /* The data is not "owned" by this structure.  It is the user's
     responsibility to free it. */
  float         *data;
};

/**
 * Initialize a lookup table to reasonable default values.
 */
int
distortion_lookup_t_init(struct distortion_lookup_t* lookup);

/**
 * Cleanup after a lookup table.  Currently does nothing, but may do
 * something in the future, so please call it when you are done with
 * the lookup table.  It does not free the data pointed to be the
 * lookup table -- it is the user's responsibility to free that array.
 */
int
distortion_lookup_t_free(struct distortion_lookup_t* lookup);

#ifdef DISTORTION_PIPELINE
/**
 * A structure to hold the parameters for a transformation based on
 * Paper IV.  Only -TAB distortions (not polynomial or spline) are
 * supported.
 *
 * Step V is performed by wcslib, so just delegates to a wcsprm
 * object.
 */
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
  double pc[NAXES*NAXES];
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
 * @return wcsini result codes
 */
int
distortion_t_init(
    struct distortion_t* dist);

/**
 * Cleanup after a distortion_t object.
 */
int
distortion_t_free(
    struct distortion_t* dist);
#endif

/**
 * Lookup the distortion offset for a particular pixel coordinate in
 * the lookup table.
 */
double
get_distortion_offset(
    const struct distortion_lookup_t *lookup,
    const double img[NAXES]);

#ifdef DISTORTION_PIPELINE
/**
 * Perform the pipeline of transformations and distortions outlined in
 * Figure 1 in Paper IV.
 *
 * Only -TAB distortions (not polynomial or spline) are
 * supported.
 *
 * @return Status codes from wcslib's wcsp2s function.
 */
int
distortion_pipeline(
    const struct distortion_t *dist,
    const unsigned int nelem,
    const double *pix /* [NAXES][nelem] */,
    /* Output parameters */
    double *world /* [NAXES][nelem] */);
#endif

/**
 * Perform just the distortion table part of Paper IV.
 *
 * @return Non-zero if an error occurred
 */
int
do_distortion(
    const unsigned int naxes,
    const struct distortion_lookup_t** lookups, /* [NAXES] */
    const unsigned int nelem,
    const double* pix, /* [NAXES][nelem] */
    double *foc /* [NAXES][nelem] */);

#endif /* __DISTORTION_H__ */
