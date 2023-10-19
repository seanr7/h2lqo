'''
Computes and plots the frequency response of a plate equipped with tuned
vibration absorbers (TVA). For more information, see
http://modelreduction.org/index.php/Plate_with_tuned_vibration_absorbers.

Copyright 2023 Quirin Aumann

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

1. Redistributions of source code must retain the above copyright notice,
this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER
OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''

from time import time
from math import pi, sqrt

import h5py
from numpy import array, vdot
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve

import matplotlib.pyplot as plt

# Specify evluation type:
#   'siso' - evaluate the z-displacement at the load location
#   'rms' - evaluate root mean square z-displacement of all nodes on the
#       plate surface
ETYPE = 'rms'

# Dimension of the matrices
n = 201900
m = 1
if ETYPE == 'siso':
    p = 1
    fname = 'plateTVA_n201900m1q1.mat'
elif ETYPE == 'rms':
    p = 28278
    n_nodes = 28278
    fname = 'plateTVA_n201900m1q28278.mat'
else:
    raise ValueError('ETYPE must be "rms" or "siso"')

# Load data
with h5py.File(fname, 'r') as f:
    def read_data(var, dtype, n1=None, n2=None):
        if dtype == 'csc_matrix':
            return csc_matrix((f.get(var + '/data')[()], \
                f.get(var + '/ir')[()], \
                f.get(var + '/jc')[()]), \
                shape=(n1, n2))
        if dtype == 'array':
            return f.get(var)[()]
        if dtype == 'cplx_array':
            return f[var]['real'] + 1j*f[var]['imag']

    B = read_data('B', 'array').T
    E = read_data('E', 'csc_matrix', n, n)
    K = read_data('K', 'csc_matrix', n, n)
    M = read_data('M', 'csc_matrix', n, n)
    s = read_data('s', 'cplx_array').flatten().tolist()

    if ETYPE == 'siso':
        C = read_data('C', 'array').T
        res = read_data('res', 'cplx_array').flatten().tolist()
    elif ETYPE == 'rms':
        C = read_data('C', 'csc_matrix', p, n)
        res = read_data('res', 'array').flatten().tolist()

# Recompute result if required
RECOMPUTE = False
if RECOMPUTE:
    res = []
    for i, mys in enumerate(s):
        print(f'Frequency step {i+1}, f={mys.imag/(2*pi):.2f} Hz ... ', end='', flush=True)

        t = time()
        A = mys**2 * M + mys * E + K
        x = spsolve(A, B)

        if ETYPE == 'siso':
            res.append(C @ x)
        elif ETYPE == 'rms':
            x = C @ x
            res.append(sqrt((vdot(x,x)).real / n_nodes))

        print(f'finished in {time() - t:.2e} s.')

# Plot the results
plt.semilogy(abs(array(s))/(2*pi), abs(array(res)))
plt.xlim([0, 250])
plt.show()

