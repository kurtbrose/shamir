'''
An implementation of shamir secret sharing algorithm.

To the extent possible under law,
all copyright and related or neighboring rights
are hereby waived.  (CC0, see LICENSE file.)

All possible patents arising from this code
are renounced under the terms of
the Open Web Foundation CLA 1.0
(http://www.openwebfoundation.org/legal/the-owf-1-0-agreements/owfa-1-0)
'''
from __future__ import division
import os
import random
import functools

# GF(256) tables will be computed at import time for the byte-based API
_GF256_EXP = [0] * 512
_GF256_LOG = [0] * 256


def _init_tables():
    '''generate log and anti-log tables for GF(256)'''
    poly = 0x11b

    def mul(a, b):
        res = 0
        for _ in range(8):
            if b & 1:
                res ^= a
            carry = a & 0x80
            a <<= 1
            if carry:
                a ^= poly
            a &= 0xff
            b >>= 1
        return res

    x = 1
    for i in range(255):
        _GF256_EXP[i] = x
        _GF256_LOG[x] = i
        x = mul(x, 3)
    for i in range(255, 512):
        _GF256_EXP[i] = _GF256_EXP[i - 255]


_init_tables()

# 12th Mersenne Prime
# (for this application we want a known prime number
# as close as possible to our security level; e.g.
# desired security level of 128 bits -- too large
# and all the ciphertext is large; too small
# and security is compromised)
_PRIME = 2**127 - 1
# 13th Mersenne Prime is 2**521 - 1

_rint = functools.partial(random.SystemRandom().randint, 0)


def _gf256_mul(a, b):
    'multiply two numbers in GF(256) using log tables'
    if a == 0 or b == 0:
        return 0
    return _GF256_EXP[_GF256_LOG[a] + _GF256_LOG[b]]


def _gf256_inv(a):
    'multiplicative inverse in GF(256)'
    if a == 0:
        raise ZeroDivisionError('inverse of 0')
    return _GF256_EXP[255 - _GF256_LOG[a]]


def _gf256_div(num, den):
    'divide in GF(256)'
    if num == 0:
        return 0
    return _gf256_mul(num, _gf256_inv(den))


def _gf256_eval_at(poly, x):
    'evaluate polynomial with coefficients in GF(256) at x'
    accum = 0
    for coeff in reversed(poly):
        accum = _gf256_mul(accum, x) ^ coeff
    return accum


def _gf256_lagrange_interpolate(x, x_s, y_s):
    'lagrange interpolation at x over GF(256)'
    k = len(x_s)
    assert k == len(set(x_s)), 'points must be distinct'
    result = 0
    for i in range(k):
        num = 1
        den = 1
        for j in range(k):
            if i == j:
                continue
            num = _gf256_mul(num, x ^ x_s[j])
            den = _gf256_mul(den, x_s[i] ^ x_s[j])
        term = _gf256_mul(y_s[i], _gf256_div(num, den))
        result ^= term
    return result


def _eval_at(poly, x, prime):
    'evaluate polynomial (coefficient tuple) at x'
    accum = 0
    for coeff in reversed(poly):
        accum *= x
        accum += coeff
        accum %= prime
    return accum


def make_random_shares(minimum, shares, prime=_PRIME):
    '''
    Generates a random shamir pool, returns
    the secret and the share points.
    '''
    if minimum > shares:
        raise ValueError("pool secret would be irrecoverable")
    poly = [_rint(prime) for i in range(minimum)]
    points = [(i, _eval_at(poly, i, prime))
              for i in range(1, shares + 1)]
    return poly[0], points


# division in integers modulus p means finding the inverse of the denominator
# modulo p and then multiplying the numerator by this inverse
# (Note: inverse of A is B such that A*B % p == 1)
# this can be computed via extended euclidean algorithm
# http://en.wikipedia.org/wiki/Modular_multiplicative_inverse#Computation
def _extended_gcd(a, b):
    x = 0
    last_x = 1
    y = 1
    last_y = 0
    while b != 0:
        quot = a // b
        a, b = b,  a%b
        x, last_x = last_x - quot * x, x
        y, last_y = last_y - quot * y, y
    return last_x, last_y


def _divmod(num, den, p):
    '''
    compute num / den modulo prime p
    To explain what this means, the return
    value will be such that the following is true:
    den * _divmod(num, den, p) % p == num
    '''
    inv, _ = _extended_gcd(den, p)
    return num * inv


def _lagrange_interpolate(x, x_s, y_s, p):
    '''
    Find the y-value for the given x, given n (x, y) points;
    k points will define a polynomial of up to kth order
    '''
    k = len(x_s)
    assert k == len(set(x_s)), "points must be distinct"
    def PI(vals):  # upper-case PI -- product of inputs
        accum = 1
        for v in vals:
            accum *= v
        return accum
    nums = []  # avoid inexact division
    dens = []
    for i in range(k):
        others = list(x_s)
        cur = others.pop(i)
        nums.append(PI(x - o for o in others))
        dens.append(PI(cur - o for o in others))
    den = PI(dens)
    num = sum([_divmod(nums[i] * den * y_s[i] % p, dens[i], p)
               for i in range(k)])
    return (_divmod(num, den, p) + p) % p


def recover_secret(shares, prime=_PRIME):
    '''
    Recover the secret from share points
    (x,y points on the polynomial)
    '''
    if len(shares) < 2:
        raise ValueError("need at least two shares")
    x_s, y_s = zip(*shares)
    return _lagrange_interpolate(0, x_s, y_s, prime)


def make_byte_shares(minimum, shares, secret):
    '''split a bytes object into share points using GF(256)'''
    if minimum > shares:
        raise ValueError('pool secret would be irrecoverable')
    if not isinstance(secret, (bytes, bytearray)):
        secret = bytes(secret)
    secret = bytearray(secret)
    polys = [
        [b] + list(os.urandom(minimum - 1))
        for b in secret
    ]
    points = []
    for x in range(1, shares + 1):
        share = bytearray(len(secret))
        for idx, poly in enumerate(polys):
            share[idx] = _gf256_eval_at(poly, x)
        points.append((x, bytes(share)))
    return points


def recover_secret_bytes(shares):
    '''recover a bytes secret from GF(256) shares'''
    if len(shares) < 2:
        raise ValueError('need at least two shares')
    x_s, y_s_bytes = zip(*shares)
    length = len(y_s_bytes[0])
    if any(len(b) != length for b in y_s_bytes):
        raise ValueError('inconsistent share lengths')
    secret = bytearray(length)
    for i in range(length):
        secret[i] = _gf256_lagrange_interpolate(0, x_s, [b[i] for b in y_s_bytes])
    return bytes(secret)


def test():
    'round trip a bunch of times; returns encrypt+decrypt time in microseconds'
    for i in range(2, 20):
        for j in range(i, i * 2):
            secret, shares = make_random_shares(i, j)
            assert recover_secret(random.sample(shares, i)) == secret
            assert recover_secret(shares) == secret
    import timeit
    return timeit.timeit(
        lambda: recover_secret(make_random_shares(4, 8)[1]),
        number=1000) * 1000


def test_gf256():
    'round trip a bunch of byte strings'
    for i in range(2, 8):
        for j in range(i, i * 2):
            secret = os.urandom(32)
            shares = make_byte_shares(i, j, secret)
            assert recover_secret_bytes(random.sample(shares, i)) == secret
            assert recover_secret_bytes(shares) == secret
    return 'ok'
