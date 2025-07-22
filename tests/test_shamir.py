import os
import random
import shamir


def test_round_trip_int():
    for i in range(2, 20):
        for j in range(i, i * 2):
            secret, shares = shamir.make_random_shares(i, j)
            assert shamir.recover_secret(random.sample(shares, i)) == secret
            assert shamir.recover_secret(shares) == secret


def test_round_trip_bytes():
    for i in range(2, 8):
        for j in range(i, i * 2):
            secret = os.urandom(32)
            shares = shamir.make_byte_shares(i, j, secret)
            assert shamir.recover_secret_bytes(random.sample(shares, i)) == secret
            assert shamir.recover_secret_bytes(shares) == secret
