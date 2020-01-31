# -*- coding: utf-8 -*-
"""
Zero forcing (Z_q)

This module implements zero forcing (Z_q version) using fast bitsets
and a brute-force approach to trying various bitsets. 
"""

#######################################################################
#
# Copyright (C) 2011 Steve Butler, Jason Grout.
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see http://www.gnu.org/licenses/.
#######################################################################

###############################################################
##### beginning of include "sage/data_structures/bitset.pxi ###
###############################################################

"""
A fast bitset datatype in Cython.

Operations between bitsets are only guaranteed to work if the bitsets
have the same size, with the exception of ``bitset_realloc``.  Similarly, you
should not try to access elements of a bitset beyond the size.

AUTHORS:

- Robert Bradshaw (2008)
- Rudi Pendavingh, Stefan van Zwam (2013-06-06): added functions map, lex_cmp,
  pickle, unpickle
- Jeroen Demeyer (2014-09-05): use mpn_* functions from MPIR in the
  implementation (:trac:`13352` and :trac:`16937`)
- Simon King (2014-10-28): ``bitset_rshift`` and ``bitset_lshift`` respecting
  the size of the given bitsets (:trac:`15820`)
"""

#*****************************************************************************
#     Copyright (C) 2008 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from libc.string cimport strlen
from cysignals.memory cimport check_calloc, check_reallocarray, sig_malloc, sig_free

from sage.cpython.string cimport char_to_str, str_to_bytes, bytes_to_str
from sage.libs.gmp.mpn cimport *
from sage.data_structures.bitset cimport *
from cython.operator import preincrement as preinc

# Doctests for the functions in this file are in sage/data_structures/bitset.pyx

#############################################################################
# Creating limb patterns
#############################################################################
#
# NOTE: In all functions in this section, the index n is interpreted
# modulo GMP_LIMB_BITS, the number of bits in a limb.
#
cdef inline mp_limb_t limb_one_set_bit(mp_bitcnt_t n):
    """
    Return a limb with only bit n set.
    """
    return (<mp_limb_t>1) << (n % GMP_LIMB_BITS)

cdef inline mp_limb_t limb_one_zero_bit(mp_bitcnt_t n):
    """
    Return a limb with all bits set, except for bit n.
    """
    return ~((<mp_limb_t>1) << (n % GMP_LIMB_BITS))

cdef inline mp_limb_t limb_lower_bits_down(mp_bitcnt_t n):
    """
    Return a limb with the lower n bits set, where n is interpreted
    in [0 .. GMP_LIMB_BITS-1].
    """
    return ((<mp_limb_t>1) << (n % GMP_LIMB_BITS)) - 1

cdef inline mp_limb_t limb_lower_bits_up(mp_bitcnt_t n):
    """
    Return a limb with the lower n bits set, where n is interpreted
    in [1 .. GMP_LIMB_BITS].
    """
    return (<mp_limb_t>(-1)) >> ((<unsigned int>(-n)) % GMP_LIMB_BITS)

#############################################################################
# Bitset Initalization
#############################################################################
cdef inline bint bitset_init(bitset_t bits, mp_bitcnt_t size) except -1:
    """
    Allocate an empty bitset of size ``size``.

    Size must be at least 1.
    """
    if size <= 0:
        raise ValueError("bitset capacity must be greater than 0")

    bits.size = size
    bits.limbs = (size - 1) // (8 * sizeof(mp_limb_t)) + 1
    bits.bits = <mp_limb_t*>check_calloc(bits.limbs, sizeof(mp_limb_t))

cdef inline int bitset_realloc(bitset_t bits, mp_bitcnt_t size) except -1:
    """
    Reallocate a bitset to size ``size``. If reallocation is larger, new bitset
    does not contain any of the extra bits.
    """
    cdef mp_size_t limbs_old = bits.limbs
    cdef mp_bitcnt_t size_old = bits.size
    if size_old == size:
        return 0
    if size <= 0:
        raise ValueError("bitset capacity must be greater than 0")

    cdef mp_size_t limbs_new = (size - 1) // (8 * sizeof(mp_limb_t)) + 1
    bits.bits = <mp_limb_t*>check_reallocarray(bits.bits, limbs_new, sizeof(mp_limb_t))
    bits.size = size
    bits.limbs = limbs_new

    if bits.limbs > limbs_old:
        # Zero any extra limbs
        mpn_zero(bits.bits + limbs_old, bits.limbs - limbs_old)
    elif bits.size < size_old:
        # Zero removed bits
        bitset_fix(bits)

cdef inline void bitset_free(bitset_t bits):
    """
    Deallocate the memory in bits.
    """
    sig_free(bits.bits)

cdef inline void bitset_clear(bitset_t bits):
    """
    Remove all elements from the set.
    """
    mpn_zero(bits.bits, bits.limbs)

cdef inline void bitset_zero(bitset_t bits):
    """
    Remove all elements from the set.

    This function is the same as bitset_clear(bits).
    """
    mpn_zero(bits.bits, bits.limbs)

cdef inline void bitset_copy(bitset_t dst, bitset_t src):
    """
    Copy the bitset src over to the bitset dst, overwriting dst.

    We assume ``dst.limbs == src.limbs``.
    """
    mpn_copyi(dst.bits, src.bits, src.limbs)

cdef inline void bitset_fix(bitset_t bits):
    """
    Clear upper bits in upper limb which should be zero.
    """
    bits.bits[bits.limbs - 1] &= limb_lower_bits_up(bits.size)

#############################################################################
# Bitset Comparison
#############################################################################

cdef inline bint mpn_equal_bits(mp_srcptr b1, mp_srcptr b2, mp_bitcnt_t n):
    """
    Return ``True`` iff the first n bits of *b1 and *b2 agree.
    """
    cdef mp_size_t nlimbs = n // GMP_LIMB_BITS
    cdef mp_limb_t mask = limb_lower_bits_down(n)
    if nlimbs > 0 and mpn_cmp(b1, b2, nlimbs) != 0:
        return False
    if mask == 0:
        return True

    cdef mp_limb_t b1h = b1[nlimbs]
    cdef mp_limb_t b2h = b2[nlimbs]
    return (b1h ^ b2h) & mask == 0

cdef bint mpn_equal_bits_shifted(mp_srcptr b1, mp_srcptr b2, mp_bitcnt_t n, mp_bitcnt_t offset):
    """
    Return ``True`` iff the first n bits of *b1 and the bits ranging from
    offset to offset+n of *b2 agree.
    """
    cdef mp_bitcnt_t bit_offset = offset % GMP_LIMB_BITS
    cdef mp_size_t i2 = offset//GMP_LIMB_BITS
    if bit_offset==0:
        return mpn_equal_bits(b1, b2 + i2, n)
    cdef mp_size_t neg_bit_offset = GMP_LIMB_BITS-bit_offset
    # limbs of b1 to be considered
    cdef mp_size_t nlimbs = n // GMP_LIMB_BITS
    # bits of an additional limb of b1 to be considered
    cdef mp_limb_t tmp_limb
    cdef mp_size_t i1
    for i1 from 0 <= i1 < nlimbs:
        tmp_limb = (b2[i2] >> bit_offset)
        tmp_limb |= (b2[preinc(i2)] << neg_bit_offset)
        if tmp_limb != b1[i1]:
            return False
    cdef mp_limb_t mask = limb_lower_bits_down(n)
    if mask == 0:
        return True

    cdef mp_limb_t b1h = b1[nlimbs]
    tmp_limb = (b2[i2] >> bit_offset)
    if (n%GMP_LIMB_BITS)+bit_offset > GMP_LIMB_BITS:
        # Need bits from the next limb of b2
        tmp_limb |= (b2[preinc(i2)] << neg_bit_offset)
    return (b1h ^ tmp_limb) & mask == 0

cdef inline bint bitset_isempty(bitset_t bits):
    """
    Test whether bits is empty.  Return True (i.e., 1) if the set is
    empty, False (i.e., 0) otherwise.
    """
    # First check lowest limb
    if bits.bits[0]:
        return False
    if bits.limbs == 1:
        return True
    # Compare bits to itself shifted by 1 limb. If these compare equal,
    # all limbs must be 0.
    return mpn_cmp(bits.bits+1, bits.bits, bits.limbs-1) == 0

cdef inline bint bitset_is_zero(bitset_t bits):
    """
    Test whether bits is empty (i.e., zero).  Return True (1) if
    the set is empty, False (0) otherwise.

    This function is the same as bitset_is_empty(bits).
    """
    return bitset_isempty(bits)

cdef inline bint bitset_eq(bitset_t a, bitset_t b):
    """
    Compare bitset a and b.  Return True (i.e., 1) if the sets are
    equal, and False (i.e., 0) otherwise.

    We assume ``a.limbs >= b.limbs``.
    """
    return mpn_cmp(a.bits, b.bits, b.limbs) == 0

cdef inline int bitset_cmp(bitset_t a, bitset_t b):
    """
    Compare bitsets a and b.  Return 0 if the two sets are
    identical, and consistently return -1 or 1 for two sets that are
    not equal.

    We assume ``a.limbs >= b.limbs``.
    """
    return mpn_cmp(a.bits, b.bits, b.limbs)

cdef inline int bitset_lex_cmp(bitset_t a, bitset_t b):
    """
    Compare bitsets ``a`` and ``b`` using lexicographical ordering.

    In this order `a < b` if, for some `k`, the first `k` elements from
    `[0 ... n-1]` are in `a` if and only if they are in `b`, and the
    `(k+1)`st element is in `b` but not ``a``. So `1010 < 1011` and
    `1010111 < 1011000`.

    We assume ``a.limbs == b.limbs``.

    INPUT:

    - ``a`` -- a bitset
    - ``b`` -- a bitset, assumed to have the same size as ``a``.

    OUTPUT:

    Return ``0`` if the two sets are identical, return ``1`` if ``a > b``,
    and return ``-1`` if ``a < b``.
    """
    cdef long i = bitset_first_diff(a, b)
    if i == -1:
        return 0
    if bitset_in(a, i):
        return 1
    else:
        return -1

cdef inline bint bitset_issubset(bitset_t a, bitset_t b):
    """
    Test whether a is a subset of b (i.e., every element in a is also
    in b).

    We assume ``a.limbs <= b.limbs``.
    """
    cdef mp_size_t i
    for i from 0 <= i < a.limbs:
        if (a.bits[i] & ~b.bits[i]) != 0:
            return False
    return True

cdef inline bint bitset_issuperset(bitset_t a, bitset_t b):
    """
    Test whether a is a superset of b (i.e., every element in b is also
    in a).

    We assume ``a.limbs >= b.limbs``.
    """
    return bitset_issubset(b, a)


cdef inline bint bitset_are_disjoint(bitset_t a, bitset_t b):
    """
    Tests whether ``a`` and ``b`` have an empty intersection.

    We assume ``a.limbs <= b.limbs``.
    """
    cdef mp_size_t i
    for i from 0 <= i < a.limbs:
        if (a.bits[i]&b.bits[i]) != 0:
            return False
    return True


#############################################################################
# Bitset Bit Manipulation
#############################################################################

cdef inline bint bitset_in(bitset_t bits, mp_bitcnt_t n):
    """
    Check if n is in bits.  Return True (i.e., 1) if n is in the
    set, False (i.e., 0) otherwise.
    """
    return (bits.bits[n >> index_shift] >> (n % GMP_LIMB_BITS)) & 1

cdef inline bint bitset_check(bitset_t bits, mp_bitcnt_t n):
    """
    Check if n is in bits.  Return True (i.e., 1) if n is in the
    set, False (i.e., 0) otherwise.

    This function is the same as bitset_in(bits, n).
    """
    return bitset_in(bits, n)

cdef inline bint bitset_not_in(bitset_t bits, mp_bitcnt_t n):
    """
    Check if n is not in bits.  Return True (i.e., 1) if n is not in the
    set, False (i.e., 0) otherwise.
    """
    return not bitset_in(bits, n)

cdef inline bint bitset_remove(bitset_t bits, mp_bitcnt_t n) except -1:
    """
    Remove n from bits.  Raise KeyError if n is not contained in bits.
    """
    if not bitset_in(bits, n):
        raise KeyError(n)
    bitset_discard(bits, n)

cdef inline void bitset_discard(bitset_t bits, mp_bitcnt_t n):
    """
    Remove n from bits.
    """
    bits.bits[n >> index_shift] &= limb_one_zero_bit(n)

cdef inline void bitset_unset(bitset_t bits, mp_bitcnt_t n):
    """
    Remove n from bits.

    This function is the same as bitset_discard(bits, n).
    """
    bitset_discard(bits, n)


cdef inline void bitset_add(bitset_t bits, mp_bitcnt_t n):
    """
    Add n to bits.
    """
    bits.bits[n >> index_shift] |= limb_one_set_bit(n)

cdef inline void bitset_set(bitset_t bits, mp_bitcnt_t n):
    """
    Add n to bits.

    This function is the same as bitset_add(bits, n).
    """
    bitset_add(bits, n)

cdef inline void bitset_set_to(bitset_t bits, mp_bitcnt_t n, bint b):
    """
    If b is True, add n to bits.  If b is False, remove n from bits.
    """
    bitset_unset(bits, n)
    bits.bits[n >> index_shift] |= (<mp_limb_t>b) << (n % GMP_LIMB_BITS)

cdef inline void bitset_flip(bitset_t bits, mp_bitcnt_t n):
    """
    If n is in bits, remove n from bits.  If n is not in bits, add n
    to bits.
    """
    bits.bits[n >> index_shift] ^= limb_one_set_bit(n)

cdef inline void bitset_set_first_n(bitset_t bits, mp_bitcnt_t n):
    """
    Set exactly the first n bits.
    """
    cdef mp_size_t i
    cdef mp_size_t index = n >> index_shift
    for i from 0 <= i < index:
        bits.bits[i] = -1
    if index < bits.limbs:
        bits.bits[index] = limb_lower_bits_down(n)
    for i from index < i < bits.limbs:
        bits.bits[i] = 0

#############################################################################
# Bitset Searching
#############################################################################

cdef inline long _bitset_first_in_limb_nonzero(mp_limb_t limb):
    """
    Given a non-zero limb of a bitset, return the index of the first
    nonzero bit.
    """
    return mpn_scan1(&limb, 0)

cdef inline long _bitset_first_in_limb(mp_limb_t limb):
    """
    Given a limb of a bitset, return the index of the first nonzero
    bit. If there are no bits set in the limb, return -1.
    """
    if limb == 0:
        return -1
    return mpn_scan1(&limb, 0)

cdef inline long bitset_first(bitset_t a):
    """
    Calculate the index of the first element in the set. If the set
    is empty, returns -1.
    """
    cdef mp_size_t i
    for i from 0 <= i < a.limbs:
        if a.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb_nonzero(a.bits[i])
    return -1

cdef inline long bitset_first_in_complement(bitset_t a):
    """
    Calculate the index of the first element not in the set. If the set
    is full, returns -1.
    """
    cdef mp_size_t i
    cdef mp_bitcnt_t j
    for i from 0 <= i < a.limbs:
        if ~a.bits[i]:
            j = (i << index_shift) | _bitset_first_in_limb_nonzero(~a.bits[i])
            if j >= a.size:
                return -1
            return <mp_size_t>j
    return -1

cdef inline long bitset_pop(bitset_t a) except -1:
    """
    Remove and return an arbitrary element from the set. Raise
    KeyError if the set is empty.
    """
    cdef long i = bitset_first(a)
    if i == -1:
        raise KeyError('pop from an empty set')
    bitset_discard(a, i)
    return i

cdef inline long bitset_first_diff(bitset_t a, bitset_t b):
    """
    Calculate the index of the first difference between a and b.  If a
    and b are equal, then return -1.

    We assume ``a.limbs == b.limbs``.
    """
    cdef mp_size_t i
    for i from 0 <= i < a.limbs:
        if a.bits[i] != b.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb_nonzero(a.bits[i] ^ b.bits[i])
    return -1

cdef inline long bitset_next(bitset_t a, mp_bitcnt_t n):
    """
    Calculate the index of the next element in the set, starting at
    (and including) n.  Return -1 if there are no elements from n
    onwards.
    """
    if n >= a.size:
        return -1
    cdef mp_size_t i = n >> index_shift
    cdef mp_limb_t limb = a.bits[i] & ~limb_lower_bits_down(n)
    cdef long ret = _bitset_first_in_limb(limb)
    if ret != -1:
        return (i << index_shift) | ret
    for i from (n >> index_shift) < i < a.limbs:
        if a.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb_nonzero(a.bits[i])
    return -1

cdef inline long bitset_next_diff(bitset_t a, bitset_t b, mp_bitcnt_t n):
    """
    Calculate the index of the next element that differs between a and
    b, starting at (and including) n.  Return -1 if there are no
    elements differing between a and b from n onwards.

    We assume ``a.limbs == b.limbs``.
    """
    if n >= a.size:
        return -1
    cdef mp_size_t i = n >> index_shift
    cdef mp_limb_t limb = (a.bits[i] ^ b.bits[i]) & ~limb_lower_bits_down(n)
    cdef long ret = _bitset_first_in_limb(limb)
    if ret != -1:
        return (i << index_shift) | ret
    for i from (n >> index_shift) < i < a.limbs:
        if a.bits[i] != b.bits[i]:
            return (i << index_shift) | _bitset_first_in_limb(a.bits[i] ^ b.bits[i])
    return -1

cdef inline long bitset_len(bitset_t bits):
    """
    Calculate the number of items in the set (i.e., the number of nonzero bits).
    """
    return mpn_popcount(bits.bits, bits.limbs)

cdef inline long bitset_hash(bitset_t bits):
    """
    Calculate a (very naive) hash function.

    This function should not depend on the size of the bitset, only on
    the items in the bitset.
    """
    cdef mp_limb_t hash = 0
    cdef mp_size_t i
    for i from 0 <= i < bits.limbs:
        hash += bits.bits[i]
    return hash

#############################################################################
# Bitset Arithmetic
#############################################################################

cdef inline void bitset_complement(bitset_t r, bitset_t a):
    """
    Set r to be the complement of a, overwriting r.

    We assume ``r.limbs == a.limbs``.
    """
    mpn_com(r.bits, a.bits, a.limbs)
    bitset_fix(r)

cdef inline void bitset_not(bitset_t r, bitset_t a):
    """
    Set r to be the complement of a, overwriting r.

    We assume ``r.limbs == a.limbs``.

    This function is the same as bitset_complement(r, a).
    """
    bitset_complement(r, a)

cdef inline void bitset_intersection(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the intersection of a and b, overwriting r.

    We assume ``a.limbs >= r.limbs == b.limbs``.
    """
    mpn_and_n(r.bits, a.bits, b.bits, b.limbs)

cdef inline void bitset_and(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the intersection of a and b, overwriting r.

    We assume ``a.limbs >= r.limbs == b.limbs``.

    This function is the same as bitset_intersection(r, a, b).
    """
    mpn_and_n(r.bits, a.bits, b.bits, b.limbs)

cdef inline void bitset_union(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the union of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    mpn_ior_n(r.bits, a.bits, b.bits, b.limbs)

cdef inline void bitset_or(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the union of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.

    This function is the same as bitset_union(r, a, b).
    """
    mpn_ior_n(r.bits, a.bits, b.bits, b.limbs)

cdef inline void bitset_difference(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the difference of a and b (i.e., things in a that are not
    in b), overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    mpn_andn_n(r.bits, a.bits, b.bits, b.limbs)

cdef inline void bitset_symmetric_difference(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the symmetric difference of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.
    """
    mpn_xor_n(r.bits, a.bits, b.bits, b.limbs)

cdef inline void bitset_xor(bitset_t r, bitset_t a, bitset_t b):
    """
    Set r to the symmetric difference of a and b, overwriting r.

    We assume ``r.limbs >= a.limbs >= b.limbs`` and either ``r is a``
    or ``r.limbs == b.limbs``.

    This function is the same as bitset_symmetric_difference(r, a, b).
    """
    mpn_xor_n(r.bits, a.bits, b.bits, b.limbs)


cdef void bitset_rshift(bitset_t r, bitset_t a, mp_bitcnt_t n):
    """
    Shift the bitset ``a`` right by ``n`` bits and store the result in
    ``r``.

    There are no assumptions on the sizes of ``a`` and ``r``.  Bits which are
    shifted outside of the resulting bitset are discarded.
    """
    if n >= a.size:
        mpn_zero(r.bits, r.limbs)
        return

    # Number of limbs on the right of a which will totally be shifted out
    cdef mp_size_t nlimbs = n >> index_shift
    # Number of limbs to be shifted assuming r is large enough
    cdef mp_size_t shifted_limbs = a.limbs - nlimbs
    # Number of bits to shift additionally
    cdef mp_bitcnt_t nbits = n % GMP_LIMB_BITS

    if shifted_limbs < r.limbs:
        if nbits:
            mpn_rshift(r.bits, a.bits + nlimbs, shifted_limbs, nbits)
        else:
            mpn_copyi(r.bits, a.bits + nlimbs, shifted_limbs)

        # Clear top limbs (note that r.limbs - shifted_limbs >= 1)
        mpn_zero(r.bits + (r.limbs - nlimbs), r.limbs - shifted_limbs)
    else:
        # Number of limbs to shift is r.limbs
        if nbits:
            mpn_rshift(r.bits, a.bits + nlimbs, r.limbs, nbits)
            if shifted_limbs > r.limbs:
                # Add the additional bits from top limb of a
                r.bits[r.limbs-1] |= a.bits[r.limbs+nlimbs] << (GMP_LIMB_BITS - nbits)
        else:
            mpn_copyi(r.bits, a.bits + nlimbs, r.limbs)

        # Clear bits outside bitset in top limb
        bitset_fix(r)

cdef void bitset_lshift(bitset_t r, bitset_t a, mp_bitcnt_t n):
    """
    Shift the bitset ``a`` left by ``n`` bits and store the result in
    ``r``.

    There are no assumptions on the sizes of ``a`` and ``r``.  Bits which are
    shifted outside of the resulting bitset are discarded.
    """
    if n >= r.size:
        mpn_zero(r.bits, r.limbs)
        return

    # Number of limbs on the right of r which will totally be zeroed
    cdef mp_size_t nlimbs = n >> index_shift
    # Number of limbs to be shifted assuming a is large enough
    cdef mp_size_t shifted_limbs = r.limbs - nlimbs
    # Number of bits to shift additionally
    cdef mp_bitcnt_t nbits = n % GMP_LIMB_BITS

    cdef mp_limb_t out = 0
    if shifted_limbs > a.limbs:
        if nbits:
            out = mpn_lshift(r.bits + nlimbs, a.bits, a.limbs, nbits)
        else:
            mpn_copyd(r.bits + nlimbs, a.bits, a.limbs)

        # Clear top limbs (note that shifted_limbs - a.limbs >= 1)
        mpn_zero(r.bits + a.limbs + nlimbs, shifted_limbs - a.limbs)
        # Store extra limb shifted in from a
        r.bits[nlimbs+a.limbs] = out
    else:
        if nbits:
            mpn_lshift(r.bits + nlimbs, a.bits, shifted_limbs, nbits)
        else:
            mpn_copyd(r.bits + nlimbs, a.bits, shifted_limbs)

        # Clear bits outside bitset in top limb
        bitset_fix(r)

    # Clear bottom limbs
    mpn_zero(r.bits, nlimbs)


cdef int bitset_map(bitset_t r, bitset_t a, m) except -1:
    """
    Fill bitset ``r`` so ``r == {m[i] for i in a}``.

    We assume ``m`` has a dictionary interface such that
    ``m[i]`` is an integer in ``[0 ... n-1]`` for all ``i`` in ``a``,
    where ``n`` is the capacity of ``r``.
    """
    cdef long i
    bitset_clear(r)
    i = bitset_first(a)
    while i >= 0:
        bitset_add(r, m[i])
        i = bitset_next(a, i + 1)
    return 0

#############################################################################
# Hamming Weights
#############################################################################

cdef inline long bitset_hamming_weight(bitset_t a):
    return bitset_len(a)

#############################################################################
# Bitset Conversion
#############################################################################

cdef char* bitset_chars(char* s, bitset_t bits, char zero=c'0', char one=c'1'):
    """
    Return a string representation of the bitset in s, using zero for
    the character representing the items not in the bitset and one for
    the character representing the items in the bitset.

    The string is both stored in s and returned.  If s is NULL, then a
    new string is allocated.
    """
    cdef mp_bitcnt_t i
    if s == NULL:
        s = <char *>sig_malloc(bits.size + 1)
    for i from 0 <= i < bits.size:
        s[i] = one if bitset_in(bits, i) else zero
    s[bits.size] = 0
    return s


cdef int bitset_from_char(bitset_t bits, char* s, char zero=c'0', char one=c'1') except -1:
    """
    Initialize a bitset with a set derived from the C string s, where one
    represents the character indicating set membership.
    """
    bitset_init(bits, strlen(s))
    cdef mp_bitcnt_t i
    for i from 0 <= i < bits.size:
        bitset_set_to(bits, i, s[i] == one)
    return 0


cdef int bitset_from_str(bitset_t bits, object s, char zero=c'0', char one=c'1') except -1:
    """
    Initialize a bitset with a set derived from the Python str s, where one
    represents the character indicating set membership.
    """
    cdef bytes b = str_to_bytes(s)
    return bitset_from_char(bits, b, zero, one)


cdef bitset_string(bitset_t bits):
    """
    Return a python string representing the bitset.
    """
    return bytes_to_str(bitset_bytes(bits))


cdef bitset_bytes(bitset_t bits):
    """
    Return a python bytes string representing the bitset.

    On Python 2 this is equivalent to bitset_string.
    """

    cdef char* s = bitset_chars(NULL, bits)
    cdef object py_s
    py_s = s
    sig_free(s)
    return py_s


cdef list bitset_list(bitset_t bits):
    """
    Return a list of elements in the bitset.
    """
    cdef list elts = []
    cdef long elt = bitset_first(bits)
    while elt >= 0:
        elts.append(elt)
        elt = bitset_next(bits, elt + 1)
    return elts

cdef bitset_pickle(bitset_t bs):
    """
    Convert ``bs`` to a reasonably compact Python structure.

    Useful for pickling objects using bitsets as internal data structure.
    To ensure this works on 32-bit and 64-bit machines, the size of a long
    is stored too.
    """
    version = 0
    data = []
    for i from 0 <= i < bs.limbs:
        data.append(bs.bits[i])
    return (version, bs.size, bs.limbs, sizeof(unsigned long), tuple(data))

cdef bitset_unpickle(bitset_t bs, tuple input):
    """
    Convert the data into a bitset.

    Companion of ``bitset_pickle()``. Assumption: ``bs`` has been initialized.
    """
    version, size, limbs, longsize, data = input
    if version != 0:
        raise TypeError("bitset was saved with newer version of Sage. Please upgrade.")
    if bs.size != size:
        bitset_realloc(bs, size)
    if sizeof(unsigned long) == longsize and bs.limbs == limbs:
        for i from 0 <= i < bs.limbs:
            bs.bits[i] = data[i]
    else:
        storage = 8 * longsize  # number of elements encoded in one limb
        adder = 0
        bitset_clear(bs)
        for i from 0 <= i < limbs:
            for j from 0 <= j < storage:
                if (data[i] >> j) & 1:
                    bitset_add(bs, j + adder)
            adder += storage

#########################################################
##### end of include "sage/data_structures/bitset.pxi ###
#########################################################

#include "sage/misc/bitset.pxi"
#from sage.misc.bitset cimport FrozenBitset, Bitset    
#include "sage/data_structures/bitset.pxi"
from sage.data_structures.bitset cimport FrozenBitset, Bitset 

### stdsage.pxi deprecated
#include "sage/ext/stdsage.pxi" 
#directly include the content intead
### Partial contents of stdsage.pxi start here ###
include "cysignals/memory.pxi"

from cysignals.memory cimport sig_malloc as sage_malloc
from cysignals.memory cimport sig_realloc as sage_realloc
from cysignals.memory cimport sig_calloc as sage_calloc
from cysignals.memory cimport sig_free as sage_free
### Partial contents of stdsage.pxi end here ###

cpdef push_zeros(list neighbors, FrozenBitset subgraph, FrozenBitset filled_set, bint return_bitset=True):
    """
    Run zero forcing as much as possible
    
    :param neighbors: (list of FrozenBitsets) -- the neighbors of each vertex
    :param subgraph: (FrozenBitset) -- the subgraph we are forcing on
    :param filled_set: (FrozenBitset) -- the initial filled vertices
    :param return_bitset: (bool) -- if True, return the set of filled vertices after playing the game,
        if False, return a boolean can_push, where can_push is True iff a force could happen.
        
    :returns: The returned filled set contains all of the initially filled vertices, even if they are outside
        of the subgraph.
    """
    cdef bitset_s *filled = &filled_set._bitset[0]
    cdef bitset_s *subgraph_bitset = &subgraph._bitset[0]

    cdef FrozenBitset v
    cdef bitset_s *current_neighbors
    cdef bitset_t unfilled_neighbors
    cdef bitset_t filled_active
    cdef bitset_t filled_active_copy
    cdef bitset_t unfilled # unfilled in the subgraph
    
    # Initialize filled_active to the complement of unfilled
    bitset_init(filled_active, filled.size)
    bitset_copy(filled_active, filled)
    
    bitset_init(unfilled, filled.size)
    bitset_complement(unfilled, filled)
    bitset_intersection(unfilled, unfilled, subgraph_bitset)

    bitset_init(filled_active_copy, filled.size)
    bitset_init(unfilled_neighbors, filled.size)
    
    cdef FrozenBitset ret = FrozenBitset(None, capacity=filled.size)
    cdef bitset_s *ret_bitset=&ret._bitset[0]

    cdef int new_filled, n
    cdef bint done = False
    cdef bint can_push = False
    
    while not done:
        done = True
        bitset_copy(filled_active_copy, filled_active)
        n = bitset_first(filled_active_copy)
        while n>=0:
            v = neighbors[n]
            current_neighbors = &v._bitset[0]
            bitset_intersection(unfilled_neighbors, current_neighbors, unfilled)
            new_filled = bitset_first(unfilled_neighbors)
            if new_filled < 0:
                # no unfilled neighbors
                bitset_discard(filled_active, n)
            else:
                # look for second unfilled neighbor
                if bitset_next(unfilled_neighbors, new_filled+1) < 0:
                    # No more unfilled neighbors
                    # push to the new_filled vertex
                    bitset_add(filled_active, new_filled)
                    bitset_remove(unfilled, new_filled)
                    bitset_remove(filled_active, n)
                    if return_bitset:
                        done = False
                    else:
                        can_push=True
                        # done is True, so the break leads to breaking out of both while loops
                        break
            n = bitset_next(filled_active_copy, n+1)

    if return_bitset:
        bitset_complement(ret_bitset, unfilled)
        bitset_intersection(ret_bitset, ret_bitset, subgraph_bitset)
        bitset_union(ret_bitset, ret_bitset, filled)

    # Free all memory used:
    bitset_free(filled_active)
    bitset_free(filled_active_copy)
    bitset_free(unfilled_neighbors)
    bitset_free(unfilled)
    
    if return_bitset:
        return ret
    else:
        return can_push
    

cpdef push_zeros_looped(list neighbors, FrozenBitset subgraph, FrozenBitset filled_set,  FrozenBitset looped, FrozenBitset unlooped, bint return_bitset=True):
    """
    Run loop zero forcing as much as possible.  Vertices that are not in the looped or unlooped sets are undetermined, so we should not apply the extra loop rules to those vertices.
    
    :param neighbors: (list of FrozenBitsets) -- the neighbors of each vertex
    :param subgraph: (FrozenBitset) -- the subgraph we are forcing on
    :param filled_set: (FrozenBitset) -- the initial filled vertices
    :param looped: (FrozenBitset) -- the vertices that are looped
    :param unlooped: (FrozenBitset) -- the vertices that are not looped
    :param return_bitset: (bool) -- if True, return the set of filled vertices after playing the game,
        if False, return a boolean can_push, where can_push is True iff a force could happen.
        
        The returned filled set contains all of the initially filled vertices, even if they are outside
        of the subgraph.
    """
    cdef bitset_s *filled = &filled_set._bitset[0]
    cdef bitset_s *subgraph_bitset = &subgraph._bitset[0]
    cdef bitset_s *unlooped_bitset = &unlooped._bitset[0]
    cdef bitset_s *looped_bitset = &looped._bitset[0]

    cdef FrozenBitset v
    cdef bitset_s *current_neighbors
    cdef bitset_t unfilled_neighbors
    cdef bitset_t active
    cdef bitset_t active_copy
    cdef bitset_t unfilled # unfilled in the subgraph
    cdef bitset_t can_die_alone
    
    
    # Initialize active to the complement of unfilled
    bitset_init(active, filled.size)
    bitset_copy(active, filled)
    # Since unlooped vertices can push early, they are "active"
    bitset_union(active, active, unlooped_bitset)
    
    bitset_init(unfilled, filled.size)
    bitset_complement(unfilled, filled)
    bitset_intersection(unfilled, unfilled, subgraph_bitset)

    bitset_init(active_copy, filled.size)
    bitset_init(unfilled_neighbors, filled.size)
    bitset_init(can_die_alone, filled.size)
    
    cdef FrozenBitset ret = FrozenBitset(None, capacity=filled.size)
    cdef bitset_s *ret_bitset=&ret._bitset[0]

    cdef int first_unfilled_neighbor, n
    cdef bint done = False
    cdef bint can_push = False
    
    while not done:
        done = True
        bitset_copy(active_copy, active)
        n = bitset_first(active_copy)
        while n>=0:
            v = neighbors[n]
            current_neighbors = &v._bitset[0]
            bitset_intersection(unfilled_neighbors, current_neighbors, unfilled)
            first_unfilled_neighbor = bitset_first(unfilled_neighbors)
            if first_unfilled_neighbor < 0:
                # no unfilled neighbors, and active means either filled or looped, so 
                # we can't die alone
                bitset_discard(active, n)
            else:
                # look for second unfilled neighbor
                if bitset_next(unfilled_neighbors, first_unfilled_neighbor+1) < 0:
                    # No more unfilled neighbors
                    # push to the first_unfilled_neighbor vertex
                    bitset_add(active, first_unfilled_neighbor)
                    bitset_remove(unfilled, first_unfilled_neighbor)
                    bitset_remove(active, n)
                    if return_bitset:
                        done = False
                    else:
                        can_push=True
                        # done is True, so the break leads to breaking out of both while loops
                        break
            n = bitset_next(active_copy, n+1)
        else:
            # only executed if we didn't break from the while n>=0 loop
            # Check to see if any looped vertex can die alone
            bitset_intersection(can_die_alone, unfilled, looped_bitset)
            n = bitset_first(can_die_alone)
            while n>=0:
                v = neighbors[n]
                current_neighbors = &v._bitset[0]
                bitset_intersection(unfilled_neighbors, current_neighbors, unfilled)
                first_unfilled_neighbor = bitset_first(unfilled_neighbors)
                if first_unfilled_neighbor < 0:
                    # no unfilled neighbors, so n dies alone
                    bitset_remove(unfilled, n)
                    if return_bitset:
                        done = False
                    else:
                        can_push=True
                        # done is True, so the break leads to breaking out of both while loops
                        break
                n = bitset_next(can_die_alone, n+1)
            
    if return_bitset:
        bitset_complement(ret_bitset, unfilled)
        bitset_intersection(ret_bitset, ret_bitset, subgraph_bitset)
        bitset_union(ret_bitset, ret_bitset, filled)

    # Free all memory used:
    bitset_free(active)
    bitset_free(active_copy)
    bitset_free(unfilled_neighbors)
    bitset_free(unfilled)
    bitset_free(can_die_alone)
    
    if return_bitset:
        return ret
    else:
        return can_push



cpdef neighbors_connected_components(list neighbors, FrozenBitset subgraph):
    cdef int n=len(neighbors)
    cdef bitset_s *subneighbors = <bitset_s *> sage_malloc(n*sizeof(bitset_s))
    cdef int i, visit
    cdef bitset_s *b
    cdef FrozenBitset FB
    cdef bitset_s *subgraph_set = &subgraph._bitset[0]
    for i in range(n):
        FB=neighbors[i]
        b=&FB._bitset[0]
        bitset_init(&subneighbors[i], n)
        bitset_intersection(&subneighbors[i], b, subgraph_set)
            
    components=set()
    cdef bitset_t seen, queue, component
    bitset_init(seen, n)
    bitset_init(queue, n)
    bitset_init(component, n)

    for i in range(n):
        if bitset_in(subgraph_set,i) and bitset_not_in(seen, i):
            bitset_clear(queue)
            bitset_clear(component)
            
            # do breadth/depth-first search from i
            bitset_add(queue,i)
            visit=i
            while visit>=0:
                bitset_remove(queue, visit)
                if bitset_not_in(component, visit):
                    bitset_add(component, visit)
                    bitset_union(queue, queue, &subneighbors[visit])
                visit=bitset_first(queue)

            components.add(tuple(bitset_list(component)))
            bitset_union(seen, seen, component)
            
            
    for i in range(n):
        bitset_free(&subneighbors[i])
    sage_free(subneighbors)

    bitset_free(seen)
    bitset_free(queue)
    bitset_free(component)
    return components
