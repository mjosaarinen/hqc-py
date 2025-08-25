#   hqc-py

2025-04-09  Markku-Juhani O. Saarinen mjos@iki.fi

2025-08-25  Updated to spec version 2025-08-22

This is a simple Python implementation of the HQC (Hamming-Quasi Cyclic)
KEM algorithm chosen by NIST for standardization in March 2025
(["winner" of the fourth round](https://www.nist.gov/news-events/news/2025/03/nist-selects-hqc-fifth-algorithm-post-quantum-encryption)).
HQC is the second PQC KEM to be standardized after Kyber/ML-KEM.
Its security is related to the Syndrome Decoding (SD) problem.

Test vectors and behavior of this implementation now matches with the
[2025-Aug-22 version of the specification](https://pqc-hqc.org/doc/hqc_specifications_2025_08_22.pdf)
from [pqc-hqc.org](https://pqc-hqc.org/) and the test vectors
of [HQC implementation](https://gitlab.com/pqc-hqc/hqc/) v5.0.0.

Here are the actual public key, secret key, ciphertext, and shared secret
sizes in bytes (_matches the reference implementation but for some
reason, not the ciphertext size in Table 7 in the current documentation._)
One can, of course, expand the secret key from a shorter seed.

| Variant |  PK  |  SK  |   CT  | SS |
|---------|------|------|-------|----|
|  HQC-1  | 2241 | 2321 |  4433 | 32 |
|  HQC-3  | 4514 | 4602 |  8978 | 32 |
|  HQC-5  | 7237 | 7333 | 14421 | 32 |

#   Code

This is a single-file implementation in [hqc.py](hqc.py), which also
includes a simple KAT generator (at the end of the file). The code only needs
SHAKE256 as a dependency from pycryptodome (`pip3 install pycryptodome`).


#   KATs and Some Unit Tests

There is an unit test file [unit_tests.py](unit_tests.py) that runs
KAT tests and also checks that the Reed-Solomon and Reed-Muller codes
are running fine (i.e. actually correcting errors).

The `Makefile` contains starts the check:

```
$ make
python3 unit_tests.py
[TEST] HQC-1 KAT 100 .................................................................................................... [END]
[PASS] HQC-1 KAT 100 84c3812eedbddde674e0a5370ecc9bfd0f71a0006cf7bcf2b1e2e26363d638a7
[TEST] HQC-3 KAT 100 .................................................................................................... [END]
[PASS] HQC-3 KAT 100 ba3f3d1e70fe73c666bede150ca7dbd0f332fc02959fe5178f8de8141b712b14
[TEST] HQC-5 KAT 100 .................................................................................................... [END]
[PASS] HQC-5 KAT 100 43dd50d6f91d9d85085558e66e2ec0168b403ded47c6dad43cd2acfddca2f618
RM HQC-1  p= 0.300  w= 5299  d= 0  True
RM HQC-1  p= 0.310  w= 5475  d= 0  True
(.. etc ...)
RS HQC-5  w= 29  True
RS HQC-5  w= 30  False
[PASS] HQC-5 RS dist 29 >= 29.
```

** EDUCATIONAL USE ONLY. ABSOLUTELY NO WARRANTY WHATSOEVER **

Cheers,

markku

