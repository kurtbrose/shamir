import random
import functools

# 12th Mersenne Prime
# (for this application we want a known prime number
# as close as possible to our security level; e.g.
# desired security level of 128 bits -- too large
# and all the ciphertext is large; too small
# and security is compromised)
_PRIME = 2**127 - 1
# 13th Mersenne Prime is 2**521 - 1

_rint = functools.partial(random.SystemRandom().randint, 0)


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
        quot = a / b
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


def recover_secret(shares):
    '''
    Recover the secret from share points
    (x,y points on the polynomial)
    '''
    pass
