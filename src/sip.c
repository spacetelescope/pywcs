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

void
sip_clear(sip_t* sip) {
  assert(sip);

  sip->a_order = 0;
  sip->a = NULL;
  sip->b_order = 0;
  sip->b = NULL;
  sip->ap_order = 0;
  sip->ap = NULL;
  sip->bp_order = 0;
  sip->bp = NULL;
  sip->crpix[0] = 0.0;
  sip->crpix[1] = 0.0;
  sip->scratch = NULL;
}

int
sip_init(
    sip_t* sip,
    const unsigned int a_order, const double* a,
    const unsigned int b_order, const double* b,
    const unsigned int ap_order, const double* ap,
    const unsigned int bp_order, const double* bp,
    const double crpix[2]) {
  unsigned int a_size       = 0;
  unsigned int b_size       = 0;
  unsigned int ap_size      = 0;
  unsigned int bp_size      = 0;
  unsigned int scratch_size = 0;
  int          status       = 0;

  assert(sip);
  sip_clear(sip);

  /* We we have one of A/B or AP/BP, we must have both. */
  if (((a == NULL) ^ (b == NULL)) ||
      ((ap == NULL) ^ (bp == NULL))) {
    return 6;
  }

  if (a != NULL) {
    sip->a_order = a_order;
    a_size = (a_order + 1) * (a_order + 1) * sizeof(double);
    sip->a = malloc(a_size);
    if (sip->a == NULL) {
      status = 2;
      goto exit;
    }
    memcpy(sip->a, a, a_size);
    if (a_order > scratch_size) {
      scratch_size = a_order;
    }

    sip->b_order = b_order;
    b_size = (b_order + 1) * (b_order + 1) * sizeof(double);
    sip->b = malloc(b_size);
    if (sip->b == NULL) {
      status = 2;
      goto exit;
    }
    memcpy(sip->b, b, b_size);
    if (b_order > scratch_size) {
      scratch_size = b_order;
    }
  }

  if (ap != NULL) {
    sip->ap_order = ap_order;
    ap_size = (ap_order + 1) * (ap_order + 1) * sizeof(double);
    sip->ap = malloc(ap_size);
    if (sip->ap == NULL) {
      status = 2;
      goto exit;
    }
    memcpy(sip->ap, ap, ap_size);
    if (ap_order > scratch_size) {
      scratch_size = ap_order;
    }

    sip->bp_order = bp_order;
    bp_size = (bp_order + 1) * (bp_order + 1) * sizeof(double);
    sip->bp = malloc(bp_size);
    if (sip->bp == NULL) {
      status = 2;
      goto exit;
    }
    memcpy(sip->bp, bp, bp_size);
    if (bp_order > scratch_size) {
      scratch_size = bp_order;
    }
  }

  if (scratch_size > 0) {
    scratch_size = (scratch_size + 1) * sizeof(double);
    sip->scratch = malloc(scratch_size);
    if (sip->scratch == NULL) {
      status = 2;
      goto exit;
    }
  }

  sip->crpix[0] = crpix[0];
  sip->crpix[1] = crpix[1];

 exit:
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
    const int x,
    const int y) {
  int index;
  assert(x >= 0 && x <= order);
  assert(y >= 0 && y <= order);
  index = x * (order + 1) + y;
  assert(index >= 0 && index < (order + 1) * (order + 1));

  return matrix[index];
}

static int
sip_compute(
    const unsigned int naxes,
    const unsigned int nelem,
    const unsigned int m,
    const double* a,
    const unsigned int n,
    const double* b,
    const double crpix[2],
    double* tmp,
    const double* input /* [NAXES][nelem] */,
    double* output /* [NAXES][nelem] */) {
  int i, j, k;
  double x, y, tmp_x, tmp_y;
  double sum;

  assert(a);
  assert(b);
  assert(crpix);
  assert(tmp);
  assert(input);
  assert(output);

  /* Avoid segfaults */
  if (input == NULL || output == NULL || tmp == NULL || crpix == NULL) {
    return 1;
  }

  /* If we have one, we must have both... */
  if ((a == NULL) ^ (b == NULL)) {
    return 6;
  }

  /* If no distortion, just copy values */
  if (a == NULL /* && b == NULL ... implied */) {
    for (i = 0; i < nelem; ++i) {
      output[i << 1] = input[i << 1];
      output[(i << 1) + 1] = input[(i << 1) + 1];
    }
    return 0;
  }

  for (i = 0; i < nelem; ++i) {
    x = input[i << 1];
    y = input[(i << 1) + 1];

    tmp_x = x - crpix[0];
    tmp_y = y - crpix[1];

    for (j = 0; j <= m; ++j) {
      tmp[j] = lu(m, a, m-j, j);
      for (k = j-1; k >= 0; --k) {
        tmp[j] = (tmp_y * tmp[j]) + lu(m, a, m-j, k);
      }
    }

    /* Don't know why this loop is convoluted */
    sum = tmp[0];
    for (j = m; j > 0; --j) {
      sum = tmp_x * sum + tmp[m - j + 1];
    }
    output[i << 1] = sum + x;

    for (j = 0; j <= n; ++j) {
      tmp[j] = lu(n, b, n-j, j);
      for (k = j-1; k >= 0; --k) {
        tmp[j] = tmp_y * tmp[j] + lu(n, b, n-j, k);
      }
    }

    /* Don't know why this loop is convoluted */
    sum = tmp[0];
    for (j = n; j > 0; --j) {
      sum = tmp_x * sum + tmp[n - j + 1];
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
                     (double *)sip->scratch,
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
                     (double *)sip->scratch,
                     foc, pix);
}

