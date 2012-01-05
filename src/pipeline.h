/*
Copyright (C) 2008-2012 Association of Universities for Research in Astronomy (AURA)

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

#ifndef __PIPELINE_H__
#define __PIPELINE_H__

#include "sip.h"
#include "distortion.h"
#include "wcs.h"

typedef struct {
  distortion_lookup_t*                   det2im[2];
  /*@shared@*/ /*@null@*/ sip_t*         sip;
  distortion_lookup_t*                   cpdis[2];
  /*@shared@*/ /*@null@*/ struct wcsprm* wcs;
  struct wcserr*                         err;
} pipeline_t;

/**
Initialize all the values in a pipeline_t to NULL.
*/
void
pipeline_clear(
    pipeline_t* pipeline);

/**
Set all the values of a pipeline_t.
*/
void
pipeline_init(
    pipeline_t* pipeline,
    /*@shared@*/ distortion_lookup_t** det2im /* [2] */,
    /*@shared@*/ sip_t* sip,
    /*@shared@*/ distortion_lookup_t** cpdis /* [2] */,
    /*@shared@*/ struct wcsprm* wcs);

/**
Free all the temporary buffers of a pipeline_t.  It does not free
the underlying sip_t, distortion_lookup_t or wcsprm objects.
*/
void
pipeline_free(
    pipeline_t* pipeline);

/**
Perform the entire pipeline from pixel coordinates to world
coordinates, in the following order:

    - Detector to image plane correction (optionally)

    - SIP distortion correction (optionally)

    - Paper IV distortion correction (optionally)

    - wcslib WCS transformation

@param ncoord:

@param nelem:

@param pixcrd [in]: Array of pixel coordinates.

@param world [out]: Array of sky coordinates (output).

@return: A wcslib error code.
*/
int
pipeline_all_pixel2world(
    pipeline_t* pipeline,
    const unsigned int ncoord,
    const unsigned int nelem,
    const double* const pixcrd /* [ncoord][nelem] */,
    double* world /* [ncoord][nelem] */);

/**
Perform just the distortion correction part of the pipeline from pixel
coordinates to focal plane coordinates.

    - Detector to image plane correction (optionally)

    - SIP distortion correction (optionally)

    - Paper IV distortion correction (optionally)

@param ncoord:

@param nelem:

@param pixcrd [in]: Array of pixel coordinates.

@param foc [out]: Array of focal plane coordinates.

@return: A wcslib error code.
*/
int
pipeline_pix2foc(
    pipeline_t* pipeline,
    const unsigned int ncoord,
    const unsigned int nelem,
    const double* const pixcrd /* [ncoord][nelem] */,
    double* foc /* [ncoord][nelem] */);

#endif
