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

#include "sip.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>

void sip_clear(sip_t* sip) {
  assert(sip);

  sip->a_order = 0;
  sip->a = NULL;
  sip->b_order = 0;
  sip->b = NULL;
  sip->ap_order = 0;
  sip->ap = NULL;
  sip->bp_order = 0;
  sip->bp = NULL;
}

int
sip_init(
    sip_t* sip,
    const unsigned int a_order, const double* a,
    const unsigned int b_order, const double* b,
    const unsigned int ap_order, const double* ap,
    const unsigned int bp_order, const double* bp,
    const double crpix[2]) {
  int status = 0;
  size_t a_size = 0;
  size_t b_size = 0;
  size_t ap_size = 0;
  size_t bp_size = 0;
  unsigned int scratch_sizes[4];
  size_t i;
  unsigned int scratch_size = 0;

  assert(sip);
  assert(sip->a == NULL);
  assert(sip->b == NULL);
  assert(sip->ap == NULL);
  assert(sip->bp == NULL);

  if (!((a == NULL) ^ (b == NULL)) ||
      !((ap == NULL) ^ (bp == NULL))) {
    return 6;
  }

  if (a != NULL) {
    sip->a_order = a_order;
    a_size = (a_order + 1) * (a_order + 1) * sizeof(double);
    sip->a = malloc(a_size);
    if (sip->a == NULL) {
      status = 2;
      goto sip_init_exit;
    }
    memcpy(sip->a, a, a_size);

    sip->b_order = b_order;
    b_size = (b_order + 1) * (b_order + 1) * sizeof(double);
    sip->b = malloc(b_size);
    if (sip->b == NULL) {
      status = 2;
      goto sip_init_exit;
    }
    memcpy(sip->b, b, b_size);
  }

  if (ap != NULL) {
    sip->ap_order = ap_order;
    ap_size = (ap_order + 1) * (ap_order + 1) * sizeof(double);
    sip->ap = malloc(ap_size);
    if (sip->ap == NULL) {
      status = 2;
      goto sip_init_exit;
    }
    memcpy(sip->ap, ap, ap_size);

    sip->bp_order = bp_order;
    bp_size = (bp_order + 1) * (bp_order + 1) * sizeof(double);
    sip->bp = malloc(bp_size);
    if (sip->bp == NULL) {
      status = 2;
      goto sip_init_exit;
    }
    memcpy(sip->bp, bp, bp_size);
  }

  sip->crpix[0] = crpix[0];
  sip->crpix[1] = crpix[1];

  /* Find the maximum order to create for our scratch memory space */
  scratch_sizes[0] = a_order + 1;
  scratch_sizes[1] = b_order + 1;
  scratch_sizes[2] = ap_order + 1;
  scratch_sizes[3] = bp_order + 1;
  for (i = 0; i < 4; ++i) {
    if (scratch_size > scratch_sizes[i]) {
      scratch_size = scratch_sizes[i];
    }
  }

  sip->scratch = malloc(scratch_size * sizeof(double));
  if (sip->scratch == NULL) {
    status = 2;
    goto sip_init_exit;
  }

 sip_init_exit:
  if (status != 0) {
    sip_free(sip);
  }

  return status;
}

void
sip_free(sip_t* sip) {
  free(sip->a);
  sip->a = NULL;
  free(sip->b);
  sip->b = NULL;
  free(sip->ap);
  sip->ap = NULL;
  free(sip->bp);
  sip->bp = NULL;
  free(sip->scratch);
  sip->scratch = NULL;
}

static inline double
lu(
    const unsigned int order,
    const double* matrix,
    const unsigned int x,
    const unsigned int y) {
  assert(x <= order);
  assert(y <= order);

  return matrix[y * (order + 1) + x];
}

static int
sip_compute(
    const unsigned int naxes,
    const unsigned int nelem,
    const unsigned int a_order,
    const double* a,
    const unsigned int b_order,
    const double* b,
    const double crpix[2],
    double* tmp,
    const double* input /* [NAXES][nelem] */,
    double* output /* [NAXES][nelem] */) {
  size_t i, j, k;
  double x, y, tmp_x, tmp_y;
  double sum;

  assert(a);
  assert(b);
  assert(crpix);
  assert(tmp);
  assert(input);
  assert(output);

  if (a == NULL || b == NULL) {
    return 6;
  }

  if (input == NULL || output == NULL) {
    return 1;
  }

  for (i = 0; i < nelem; ++i) {
    x = input[i << 1] + 1.0;
    y = input[(i << 1) + 1] + 1.0;

    tmp_x = crpix[0] - x;
    tmp_y = crpix[1] - y;

    for (j = 0; j <= a_order; ++j) {
      tmp[j] = lu(a_order, a, a_order - j, j);
      for (k = j - 1; k >= 0; --k) {
        tmp[j] = tmp_y * tmp[j] + lu(a_order, a, a_order - j, k);
      }
    }

    sum = tmp[0];
    for (i = a_order; i > 0; --i) {
      sum = tmp_x * sum + tmp[a_order - i + 1];
    }
    output[i << 1] = sum + x;

    for (j = 0; j <= b_order; ++j) {
      tmp[j] = lu(b_order, b, b_order - j, j);
      for (k = j - 1; k >= 0; --k) {
        tmp[j] = tmp_y * tmp[j] + lu(b_order, b, b_order - j, k);
      }
    }

    sum = tmp[0];
    for (i = b_order; i > 0; --i) {
      sum = tmp_x * sum + tmp[b_order - i + 1];
    }
    output[(i << 1) + 1] = sum + y;
  }

  return 0;
}

int
sip_pix2foc(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* pix /* [NAXES][nelem] */,
    double* foc /* [NAXES][nelem] */) {
  if (sip == NULL) {
    return 1;
  }

  return sip_compute(naxes, nelem,
                     sip->a_order, sip->a,
                     sip->b_order, sip->b,
                     sip->crpix,
                     sip->scratch,
                     pix, foc);
}

int
sip_foc2pix(
    const sip_t* sip,
    const unsigned int naxes,
    const unsigned int nelem,
    const double* foc /* [NAXES][nelem] */,
    double* pix /* [NAXES][nelem] */) {
  if (sip == NULL) {
    return 1;
  }

  return sip_compute(naxes, nelem,
                     sip->ap_order, sip->ap,
                     sip->bp_order, sip->bp,
                     sip->crpix,
                     sip->scratch,
                     foc, pix);
}

