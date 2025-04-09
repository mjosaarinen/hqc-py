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
includes a KAT generator (at the end of the file). The code only needs
SHAKE256 as a dependency from pycryptodome (`pip3 install pycryptodome`).


I've include sha256 hashes of the test vectors for HQC-128, HQC-192,
HQC-256 only as those are very large. The `Makefile` contains a simple
check:

```
$ make
grep `python3 hqc.py | sha256sum | colrm 65` kat/*
kat/hqc-kat-hash.txt:2be3afb1efb98ce58d719e19824f1a1fcb53fedb6ca8e16bc34afa53ac0528c0  hqc-128.rsp.100
```
To change the test target (HQC-192 etc), change the last line in `hqc.py`

I've also unit-tested the decoders with an artificially large number of
errors (Reed-Solomun up to "delta" errors, etc.), and those should be
working fine. The decoding algorithms are, of course, very elementary ones
in this implementation, and no effort is made for constant time, etc.

** EDUCATIONAL USE ONLY. ABSOLUTELY NO WARRANTY WHATSOEVER **

Cheers,

markku

