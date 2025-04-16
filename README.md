#   hqc-py

2025-04-09  Markku-Juhani O. Saarinen mjos@iki.fi

This is a simple Python implementation of the HQC (Hamming-Quasi Cyclic)
KEM algorithm chosen by NIST for standardization in March 2025
(["winner" of the fourth round](https://www.nist.gov/news-events/news/2025/03/nist-selects-hqc-fifth-algorithm-post-quantum-encryption)).
HQC is the second PQC KEM to be standardized after Kyber/ML-KEM.
Its security is related to the Syndrome Decoding (SD) problem.

Test vectors and behavior should match with the 2025-Feb-19 version of the specification from  [pqc-hqc.org](https://pqc-hqc.org/).

Here are the actual public key, secret key, ciphertext, and shared secret
sizes in bytes (_matches the reference implementation but for some
reason, not the ciphertext size in Table 7 in the current documentation._)
One can, of course, expand the secret key from a shorter seed.

| Variant |  PK  |  SK  |   CT  | SS |
|---------|------|------|-------|----|
| hqc-128 | 2249 | 2305 |  4433 | 64 |
| hqc-192 | 4522 | 4586 |  8978 | 64 |
| hqc-256 | 7245 | 7317 | 14421 | 64 |

#   Code

This is a single-file implementation in [hqc.py](hqc.py), which also
includes a simple KAT generator (at the end of the file). The code only needs
SHAKE256 as a dependency from pycryptodome (`pip3 install pycryptodome`).


#	KATs and Some Unit Tests

There is an unit test file [unit_tests.py](unit_tests.py) that runs
KAT tests and also checks that the Reed-Solomon and Reed-Muller codes
are running fine (i.e. actually correcting errors).

The `Makefile` contains starts the check:

```
$ make
python3 unit_tests.py
[TEST] hqc-128 KAT 100 .................................................................................................... [END]
[PASS] hqc-128 KAT 100 2be3afb1efb98ce58d719e19824f1a1fcb53fedb6ca8e16bc34afa53ac0528c0
[TEST] hqc-192 KAT 100 .................................................................................................... [END]
[PASS] hqc-192 KAT 100 a77dd028b01a6554eebe5530fbccf4891e788c50a7efae9156e7ebc01bbea20a
[TEST] hqc-256 KAT 100 .................................................................................................... [END]
[PASS] hqc-256 KAT 100 a63d941b042b3992dc9d19f9b901548bfcd7e5b9cb92bb87724ccfbbabee8917
RM hqc-128  p= 0.300  w= 5299  d= 0  True
RM hqc-128  p= 0.310  w= 5475  d= 0  True
(.. etc ...)
RS hqc-256  w= 29  True
RS hqc-256  w= 30  False
[PASS] hqc-256 RS dist 29 >= 29.
```

** EDUCATIONAL USE ONLY. ABSOLUTELY NO WARRANTY WHATSOEVER **

Cheers,

markku

