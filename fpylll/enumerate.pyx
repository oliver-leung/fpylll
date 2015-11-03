# -*- coding: utf-8 -*-
# distutils: language = c++
# distutils: libraries = gmp mpfr fplll

include "interrupt/interrupt.pxi"

from libcpp.vector cimport vector
from gso cimport MatGSO
from fplll cimport Enumeration as Enumeration_c
from fplll cimport MatGSO as MatGSO_c
from fplll cimport Z_NR, FP_NR, mpz_t
from fplll cimport FastEvaluator

cdef class Enumeration:

    @staticmethod
    def enumerate(MatGSO m, max_dist, max_dist_expo, first, last, pruning):

        cdef MatGSO_c[Z_NR[mpz_t], FP_NR[double]] *gso = m._core_mpz_double

        cdef FastEvaluator[FP_NR[double]] evaluator

        cdef vector[double] pruning_
        cdef vector[FP_NR[double]] target_coord_
        cdef vector[FP_NR[double]] sub_tree_

        cdef int block_size = last-first

        if pruning is None:
            for i in range(block_size):
                pruning_.push_back(1)
        else:
            for i in range(block_size):
                pruning_.push_back(pruning[i])

        cdef double max_dist__ = max_dist
        cdef FP_NR[double] max_dist_ = max_dist__

        sig_on()
        Enumeration_c.enumerate[FP_NR[double]](gso[0], max_dist_, max_dist_expo, evaluator,
                                               target_coord_, sub_tree_,
                                               first, last, pruning_)
        sig_off()

        solution = []
        for i in range(evaluator.solCoord.size()):
            solution.append(evaluator.solCoord[i].get_d())

        max_dist = max_dist_.get_d()
        return tuple(solution), max_dist