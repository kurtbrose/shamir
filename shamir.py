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


def _lagrange_interpolate(x, x_s, y_s):
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
    num = sum([nums[i] * den * y_s[i] / dens[i] for i in range(k)])
    assert num % den == 0, "bad points: got non-integer y value"
    return num / den


def recover_secret(shares, prime=_PRIME):
    '''
    Recover the secret from share points
    (x,y points on the polynomial)
    '''
    if len(shares) < 2:
        raise ValueError("need at least two shares")
    x_s, y_s = zip(*shares)
    return _lagrange_interpolate(0, x_s, y_s) % prime


def test():
    for i in range(3, 20):
        for j in range(i, i * 2):
            secret, shares = make_random_shares(i, j)
            assert recover_secret(random.sample(shares, i)) == secret
            assert recover_secret(shares) == secret
            print '.',
    import timeit
    print timeit.timeit(
        lambda: recover_secret(make_random_shares(4, 8)),
        number=1000) * 1000