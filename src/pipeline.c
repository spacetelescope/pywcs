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
pipeline_clear(
    pipeline_t* pipeline) {

  pipeline->sip = NULL;
  pipeline->cpdis[0] = NULL;
  pipeline->cpdis[1] = NULL;
  pipeline->wcs = NULL;

  /* Temporary buffers */
  pipeline->alloc_ncoord = 0;
  pipeline->alloc_nelem = 0;
  pipeline->tmp = NULL;
  pipeline->imgcrd = NULL;
  pipeline->phi = NULL;
  pipeline->theta = NULL;
  pipeline->stat = NULL;
}

void
pipeline_init(
    pipeline_t* pipeline,
    /*@shared@*/ sip_t* sip,
    /*@shared@*/ distortion_lookup_t** cpdis /* [2] */,
    /*@shared@*/ struct wcsprm* wcs) {

  pipeline->sip = sip;
  pipeline->cpdis[0] = cpdis[0];
  pipeline->cpdis[1] = cpdis[1];
  pipeline->wcs = wcs;
}

static void
pipeline_free_tmp(
    pipeline_t* pipeline) {

  /* Free all temporary buffers and reset pointers to NULL */
  pipeline->alloc_nelem = 0;
  pipeline->alloc_ncoord = 0;
  free(pipeline->tmp);
  pipeline->tmp = NULL;
  free(pipeline->imgcrd);
  pipeline->imgcrd = NULL;
  free(pipeline->phi);
  pipeline->phi = NULL;
  free(pipeline->theta);
  pipeline->theta = NULL;
  free(pipeline->stat);
  pipeline->stat = NULL;
}

void
pipeline_free(
    pipeline_t* pipeline) {

  pipeline_free_tmp(pipeline);
}

static int
pipeline_realloc(
    pipeline_t* pipeline,
    unsigned int ncoord,
    unsigned int nelem) {

  if (pipeline->alloc_ncoord < ncoord ||
      pipeline->alloc_nelem * pipeline->alloc_ncoord < nelem * ncoord) {
    pipeline_free_tmp(pipeline);

    pipeline->imgcrd = malloc(ncoord * nelem * sizeof(double));
    if (pipeline->imgcrd == NULL) {
      goto out_of_memory;
    }

    pipeline->phi = malloc(ncoord * sizeof(double));
    if (pipeline->phi == NULL) {
      goto out_of_memory;
    }

    pipeline->theta = malloc(ncoord * sizeof(double));
    if (pipeline->theta == NULL) {
      goto out_of_memory;
    }

    pipeline->stat = malloc(ncoord * nelem * sizeof(int));
    if (pipeline->stat == NULL) {
      goto out_of_memory;
    }

    pipeline->tmp = malloc(ncoord * nelem * sizeof(double));
    if (pipeline->tmp == NULL) {
      goto out_of_memory;
    }

    pipeline->alloc_ncoord = ncoord;
    pipeline->alloc_nelem = nelem;
  }

  return 0;

 out_of_memory:
  pipeline->alloc_nelem = 0;
  pipeline->alloc_ncoord = 0;
  pipeline_free_tmp(pipeline);
  return 2;
}

int
pipeline_all_pixel2world(
    const pipeline_t* pipeline,
    const unsigned int ncoord,
    const unsigned int nelem,
    const double* pixcrd /* [ncoord][nelem] */,
    double* world /* [ncoord][nelem] */) {

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
    status = pipeline_realloc((pipeline_t*)pipeline, ncoord, nelem);
    if (status != 0) {
      goto exit;
    }

    if (has_sip || has_p4) {
      status = pipeline_pix2foc(pipeline, ncoord, nelem, pixcrd, pipeline->tmp);
      if (status != 0) {
        goto exit;
      }

      wcs_input = pipeline->tmp;
      wcs_output = world;
    } else {
      wcs_input = pixcrd;
      wcs_output = world;
    }

    status = wcsp2s(pipeline->wcs, (int)ncoord, (int)nelem,
                    wcs_input, pipeline->imgcrd, pipeline->phi,
                    pipeline->theta, wcs_output, pipeline->stat);
  } else {
    if (has_sip || has_p4) {
      status = pipeline_pix2foc(pipeline, ncoord, nelem, pixcrd, world);
    }
  }

 exit:
  /* We don't have any cleanup at the moment */
  return status;
}

int pipeline_pix2foc(
    const pipeline_t* pipeline,
    const unsigned int ncoord,
    const unsigned int nelem,
    const double* pixcrd /* [ncoord][nelem] */,
    double* foc /* [ncoord][nelem] */) {

  int           has_sip;
  int           has_p4;
  int           status     = 0;

  assert(nelem == 2);

  if (pipeline == NULL || pixcrd == NULL || foc == NULL) {
    return 1;
  }

  has_sip = pipeline->sip != NULL;
  has_p4  = pipeline->cpdis[0] != NULL && pipeline->cpdis[1] != NULL;

  /* Copy pixcrd to foc as a starting point.  The "deltas" functions below will
     undistort from there */
  memcpy(foc, pixcrd, sizeof(double) * ncoord * nelem);

  if (has_sip) {
    status = sip_pix2deltas(pipeline->sip, 2, ncoord, pixcrd, foc);
    if (status) {
      goto exit;
    }
  }

  if (has_p4) {
    status = p4_pix2deltas(2, (void*)pipeline->cpdis, ncoord, pixcrd, foc);
    if (status) {
      goto exit;
    }
  }

 exit:
  /* We don't have any cleanup at the moment */
  return status;
}

