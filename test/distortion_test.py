import pywcs
import numpy as np

d1 = np.asarray(np.random.random((50,50)), dtype=np.float32) - 0.5
d2 = np.asarray(np.random.random((50,50)), dtype=np.float32) - 0.5
cpdis_x = pywcs.DistortionLookupTable(d1,
                                      (25, 25), (25, 25), (25, 25))
cpdis_y = pywcs.DistortionLookupTable(d2,
                                      (25, 25), (25, 25), (25, 25))
cpdis_z = pywcs.DistortionLookupTable(d1,
                                      (25, 25), (25, 25), (25, 25))
pos = np.array([[1024,1024]])
print pywcs._pywcs.do_distortion([cpdis_x, cpdis_y], pos)
print pos

print cpdis_x.data
print cpdis_y.data
del cpdis_z

# dist = pywcs.Distortion()
# dist.cd = np.array([[1, 0], [0, 1]])
# dist.crpix = [-234.75, 8.3393]
# dist.cdelt = np.array([-0.066667, 0.066667])
# dist.crval = [0, -90]
# dist.ctype = ["RA---TAN", "DEC--TAN"]
# dist.cpdis = (cpdis_x, cpdis_y)
# print dist.crval
# print dist.ctype
# print dist.p2s(np.random.random((2, 50)))
