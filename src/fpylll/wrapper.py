# -*- coding: utf-8 -*-
"""
Wrapper.

..  moduleauthor:: Oliver Leung <oliverleung2@gmail.com>

This class provides an externally callable API for LLL reducing some basis ``b``.
This class forcibly instantiates some template declarations (see FPLLL_DECLARE_LLL(T) for more
information at the bottom of this class) and provides an interface for calling provable,
heuristic and fast variants of both LLL and HLLL.

User note: all of the parameters passed into this class are non-const references. Thus, this
class may overwrite these variables.

In particular, the methods in this class typically
over-write the parameter ``u`` with the transformation matrix applied to ``b`` to produce LLL(B).
In other words, let C = LLL(B). Then we can write C = U B for some unimodular transformation
matrix U. For this reason, ``u`` should either be an empty matrix or the identity matrix when it
is passed in to this function, so that no information is lost.

Similarly, the parameter ``u_inv`` will contain the inverse matrix ``u ^ -1``.
This allows you to recover the original basis ``b`` by multiplying ``C`` by ``u_inv``.
Note that all operations on this object are carried out on the transpose of this matrix
(see gso_interface.h) for speed: hence the discrepancy between the names (see the lll_reduction
method in wraper.cpp for this). This and other such behaviour on these parameters is described in
some detail in gso_interface.h.
"""
from math import ceil, log

from fpylll import IntegerMatrix
from fpylll.util import get_precision

DIM_DOUBLE_MAX = [
    0,     26,    29.6,  28.1,  31.1,  32.6,  34.6,  34,    37.7,  38.8,  39.6,  41.8,  40.9,
    43.6,  44.2,  47,    46.8,  50.6,  49.1,  51.5,  52.5,  54.8,  54.6,  57.4,  57.6,  59.9,
    61.8,  62.3,  64.5,  67.1,  68.8,  68.3,  69.9,  73.1,  74,    76.1,  76.8,  80.9,  81.8,
    83,    85.3,  87.9,  89,    90.1,  89,    94.6,  94.8,  98.7,  99,    101.6, 104.9, 106.8,
    108.2, 107.4, 110,   112.7, 114.6, 118.1, 119.7, 121.8, 122.9, 126.6, 128.6, 129,   133.6,
    126.9, 135.9, 139.5, 135.2, 137.2, 139.3, 142.8, 142.4, 142.5, 145.4]

ETA_DEP = [1.,       # 0.5
           1.,       # 0.55
           1.0521,   # 0.6
           1.1254,   # 0.65
           1.2535,   # 0.7
           1.3957,   # 0.75
           1.6231,   # 0.8
           1.8189,   # 0.85
           2.1025,   # 0.9
           2.5117]  # 0.95

class Wrapper:

    def __init__(self, b: IntegerMatrix, u, u_inv, delta, eta, flags, theta=None, c=None):
        # self.status = RED_SUCCESS # How to import this from fplll.pxd?
        self.b = b
        self.u = u
        self.u_inv = u_inv
        self.delta = delta
        self.eta = eta
        self.use_long = False

        self.n = b.ncols
        self.d = b.nrows
        self.flags = flags

        if theta is None and c is None:
            self.last_early_red = 0
            self.max_exponent = b.get_max_exp() + ceil(0.5 * log(self.d * self.n, 2))
            self.good_prec = 10 # how to check min prec?
            # self.good_prec = get_precision()

        elif theta is not None and c is not None:
            self.last_early_red = -1
            self.good_prec = 10  # how to check min prec?
            # self.good_prec = get_precision()

    def little(self, kappa, precision):
        dm = int(self.delta * 100 - 25)
        if dm < 0:
            dm = 0
        elif dm >= len(DIM_DOUBLE_MAX):
            dm = len(DIM_DOUBLE_MAX) - 1

        em = int((self.eta - 0.5) * 20)
        if em < 0:
            em = 0
        elif em >= len(ETA_DEP):
            em = len(ETA_DEP) - 1

        p = max(1, precision / 53)
        p *= DIM_DOUBLE_MAX[dm]
        p *= ETA_DEP[em]

        return kappa < p
