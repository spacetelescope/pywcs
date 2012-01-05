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

#define NO_IMPORT_ARRAY

#include "pywcs_api.h"

int
PyWcs_GetCVersion(void) {
  return REVISION;
}

void* PyWcs_API[] = {
  /*  0 */ (void *)PyWcs_GetCVersion,
  /* pyutil.h */
  /*  1 */ (void *)wcsprm_python2c,
  /*  2 */ (void *)wcsprm_c2python,
  /* distortion.h */
  /*  3 */ (void *)distortion_lookup_t_init,
  /*  4 */ (void *)distortion_lookup_t_free,
  /*  5 */ (void *)get_distortion_offset,
  /*  6 */ (void *)p4_pix2foc,
  /*  7 */ (void *)p4_pix2deltas,
  /* sip.h */
  /*  8 */ (void *)sip_clear,
  /*  9 */ (void *)sip_init,
  /* 10 */ (void *)sip_free,
  /* 11 */ (void *)sip_pix2foc,
  /* 12 */ (void *)sip_pix2deltas,
  /* 13 */ (void *)sip_foc2pix,
  /* 14 */ (void *)sip_foc2deltas,
  /* pipeline.h */
  /* 15 */ (void *)pipeline_clear,
  /* 16 */ (void *)pipeline_init,
  /* 17 */ (void *)pipeline_free,
  /* 18 */ (void *)pipeline_all_pixel2world,
  /* 19 */ (void *)pipeline_pix2foc,
  /* wcs.h */
  /* 20 */ (void *)wcsp2s,
  /* 21 */ (void *)wcss2p,
  /* 22 */ (void *)wcsprt,
  /* new for api version 2 */
  /* 23 */ (void *)wcslib_get_error_message,
  /* new for api version 3 */
  /* 24 */ (void *)wcsprintf_buf
};

int _setup_api(PyObject *m) {
  PyObject* c_api;

  #if PY_VERSION_HEX >= 0x03020000
    c_api = PyCapsule_New((void *)PyWcs_API, "_pywcs._PYWCS_API", NULL);
  #else
    c_api = PyCObject_FromVoidPtr((void *)PyWcs_API, NULL);
  #endif
  PyModule_AddObject(m, "_PYWCS_API", c_api);

  return 0;
}
