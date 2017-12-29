shamir
''''''

Fast, secure, pure python implementation of Shamir's secret sharing algorithm.

https://en.wikipedia.org/wiki/Shamir%27s_Secret_Sharing


why to use it
'''''''''''''

Shamir's secret sharing is useful for providing an "N of M" layer.
For example, password recovery security questions could be cryptographically
imeplemented with this algorithm.  Security question answers are
put through a key-derivation-function (KDF) or hash and each one is used to
encrypt a different secret.  Then, the answer to any say 5 of 11 security
questions would be enouch to recover the secret.

how to use it
'''''''''''''

.. code-block:: python
    secret, shares = make_random_secret(3, 5)
    # generate shares such that 3 of 5 can recover the secret
    secret = recover_secret(shares)


Shamir's secret sharing algorithm operates on integer X, Y points,
and the secret it stores is a random integer.  To be useful, it must be
combined with other algorithms.  Here's a high level example:

.. code-block:: python
    # encrypt, decrypt, hash, and kdf are external functions

    def two_of_three_encrypt(plaintext, pw0, pw1, pw2):
        'given a plaintext, secure it so that any 2 may access in the future'
        secret, shares = shamir.make_random_shares(2, 3)
        def encrypt_share(share, pw):
            return encrypt(key=kdf(pw), plaintext=repr(share))
        return (
            encrypt(key=hash(hex(secret)), plaintext=plaintext)) + tuple(
            [encrypt(key=kef(pw), plaintext=repr(share))
             for pw, share in zip((pw0, pw1, pw2), shares)])

    def two_of_three_decrypt(encrypted, pwA, pwB):
        'recover the plaintext given 2 of the 3 passwords used to secure'
        ciphertext, shares = encrypted[0], encrypted[1:]
        keyA, keyB = kdf(pwA), kdf(pwB)
        decrypted_shares = []
        for share in shares:
            for key in (keyA, keyB):
                try:
                    decrypted_shares.append(
                        decrypt(key=keyA, ciphertext=share))
                except Exception:
                    pass
        if len(decrypted_shares) < 2:
            raise ValueError('bad password')
        return decrypt(
            key=hash(hex(shamir.recover_secret(decrypted_shares))),
            ciphertext=ciphertext)
