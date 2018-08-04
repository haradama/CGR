# coding: utf-8

import numpy as np
cimport numpy as np
import itertools

DTYPE_INT = np.int
DTYPE_FLOAT = np.float

ctypedef np.int_t DTYPE_INT_t
ctypedef np.float_t DTYPE_FLOAT_t

cdef class CGR:
    cdef int k, array_size
    def __cinit__(self, int k_length):        
        self.k = k_length
        self.array_size = int((4**k_length)**0.5)

    cpdef np.ndarray[DTYPE_FLOAT_t, ndim=2] cgr(self, str seq):
        cdef str nucl, fragment
        cdef int index, maxX, maxY, locX, locY, seq_length

        seq_length = len(seq)

        cdef np.ndarray[DTYPE_INT_t, ndim=2] chaos = np.zeros([self.array_size, self.array_size], dtype=DTYPE_INT)

        maxX, maxY = self.array_size, self.array_size
        locX, locY = 1, 1

        for index in range(seq_length - (self.k - 1)):
            fragment = seq[index:index + self.k]
            if "N" not in fragment:
                for nucl in fragment:
                    if nucl == "T":
                        locX += maxX // 2
                    elif nucl == "C":
                        locY += maxY // 2
                    elif nucl == "G":
                        locX += maxX // 2
                        locY += maxY // 2

                    maxX //= 2
                    maxY //= 2

                chaos[locY-1, locX-1] += 1
                maxX = self.array_size
                maxY = self.array_size
                locX = 1
                locY = 1

        return chaos / (seq_length - (self.k - 1))
    
    cpdef get_signatures(self):
        signatures = np.full((self.array_size, self.array_size), "N" * self.k)

        cdef int maxX, maxY, locX, locY

        base = ("A", "T", "G", "C")

        maxX, maxY = self.array_size, self.array_size
        locX, locY = 1, 1

        for p in itertools.product(*([base] * self.k)):
            for nucl in p:
                if nucl == "T":
                    locX += maxX // 2
                elif nucl == "C":
                    locY += maxY // 2
                elif nucl == "G":
                    locX += maxX // 2
                    locY += maxY // 2

                maxX //= 2
                maxY //= 2

            signatures[locY-1, locX-1] = "".join(p)
            maxX = self.array_size
            maxY = self.array_size
            locX = 1
            locY = 1

        return signatures

cpdef double GC(str seq):
    return (seq.count("G") + seq.count("C")) / len(seq)