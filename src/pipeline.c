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

#include "pipeline.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>

void
pipeline_clear(pipeline_t* pipeline) {
  pipeline->sip = NULL;
  pipeline->cpdis[0] = NULL;
  pipeline->cpdis[1] = NULL;
  pipeline->wcs = NULL;
}

void
pipeline_init(
    pipeline_t* pipeline,
    sip_t* sip,
    distortion_lookup_t* cpdis[2],
    struct wcsprm* wcs) {
  pipeline->sip = sip;
  pipeline->cpdis[0] = cpdis[0];
  pipeline->cpdis[1] = cpdis[1];
  pipeline->wcs = wcs;
}

int
pipeline_all_pixel2world(
    const pipeline_t* pipeline,
    const int ncoord,
    const int nelem,
    const double* pixcrd /* [ncoord][nelem] */,
    double* world /* [ncoord][nelem] */) {
  double*       tmp        = NULL;
  double*       imgcrd     = NULL;
  double*       phi        = NULL;
  double*       theta      = NULL;
  int*          stat       = NULL;
  const double* wcs_input  = NULL;
  double*       wcs_output = NULL;
  int           has_sip;
  int           has_p4;
  int           has_wcs;
  int           status     = 0;

  assert(nelem == 2);

  if (pipeline == NULL || pixcrd == NULL || world == NULL) {
    return 1;
  }

  has_sip = pipeline->sip != NULL;
  has_p4  = pipeline->cpdis[0] != NULL && pipeline->cpdis[1] != NULL;
  has_wcs = pipeline->wcs != NULL;

  if (has_wcs) {
    imgcrd = malloc(ncoord * nelem * sizeof(double));
    if (imgcrd == NULL) {
      status = 2;
      goto pipeline_pixel2world_exit;
    }

    phi = malloc(ncoord * sizeof(double));
    if (phi == NULL) {
      status = 2;
      goto pipeline_pixel2world_exit;
    }

    theta = malloc(ncoord * sizeof(double));
    if (theta == NULL) {
      status = 2;
      goto pipeline_pixel2world_exit;
    }

    stat = malloc(ncoord * nelem * sizeof(int));
    if (stat == NULL) {
      status = 2;
      goto pipeline_pixel2world_exit;
    }

    if (has_sip || has_p4) {
      tmp = malloc(ncoord * nelem * sizeof(double));
      if (tmp == NULL) {
        status = 2;
        goto pipeline_pixel2world_exit;
      }

      status = pipeline_pix2foc(pipeline, ncoord, nelem, pixcrd, tmp);
      if (status != 0) {
        goto pipeline_pixel2world_exit;
      }

      wcs_input = tmp;
      wcs_output = world;
    } else {
      wcs_input = pixcrd;
      wcs_output = world;
    }

    status = wcsp2s(pipeline->wcs, ncoord, nelem,
                    wcs_input, imgcrd, phi, theta, wcs_output, stat);
  } else {
    if (has_sip || has_p4) {
      status = pipeline_pix2foc(pipeline, ncoord, nelem, pixcrd, world);
    }
  }

 pipeline_pixel2world_exit:
  free(tmp);
  free(imgcrd);
  free(phi);
  free(theta);
  free(stat);

  return status;
}

int pipeline_pix2foc(
    const pipeline_t* pipeline,
    const int ncoord,
    const int nelem,
    const double* pixcrd /* [ncoord][nelem] */,
    double* foc /* [ncoord][nelem] */) {
  const double* sip_input  = NULL;
  double*       sip_output = NULL;
  const double* p4_input   = NULL;
  double*       p4_output  = NULL;
  int           has_sip;
  int           has_p4;
  int           status     = 0;

  assert(nelem == 2);

  if (pipeline == NULL || pixcrd == NULL || foc == NULL) {
    return 1;
  }

  has_sip = pipeline->sip != NULL;
  has_p4  = pipeline->cpdis[0] != NULL && pipeline->cpdis[1] != NULL;

  /* Build a plan for how to use the buffers. */
  if (has_sip) {
    if (has_p4) {
      sip_input  = pixcrd;
      sip_output = foc;
      p4_input   = foc;
      p4_output  = foc;
    } else {
      sip_input  = pixcrd;
      sip_output = foc;
    }
  } else {
    if (has_p4) {
      p4_input  = pixcrd;
      p4_output = foc;
    } else {
      return 0;
    }
  }

  if (has_sip) {
    status = sip_pix2foc(pipeline->sip, 2, ncoord, sip_input, sip_output);
    if (status) {
      goto pipeline_pix2foc_exit;
    }
  }

  if (has_p4) {
    status = p4_pix2foc(2, (void*)pipeline->cpdis, ncoord, p4_input, p4_output);
    if (status) {
      goto pipeline_pix2foc_exit;
    }
  }

 pipeline_pix2foc_exit:
  /* We don't have any cleanup at the moment */

  return status;
}

