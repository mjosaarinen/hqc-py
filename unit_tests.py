#   unit_tests.py
#   2025-04-03  Markku-Juhani O. Saarinen <mjos@iki.fi>

#   === Unit test for error correction

from Crypto.Hash import SHAKE256, SHA256
import random
from hqc import HQC_ALL

#   sha-256 hashes of KAT rsp files -- also limited to 1,10 first ones
HQC_KAT_SHA256 = {
  (1, 1):   'f55cf152d07ff1b53c994bfd078c634830edbbbc39d8b62ceab26c80da201148',
  (1, 10):  '8eb9432587f0efa5f01a12d8cd3be940bf1673dd37a7d0299299f0a51651c1bf',
  (1, 100): '84c3812eedbddde674e0a5370ecc9bfd0f71a0006cf7bcf2b1e2e26363d638a7',
  (3, 1):   '14ef57363010458280e9e78a9a6aada7167f32fc90398663206a9334e2d54e56',
  (3, 10):  '7d933ee869815ef5fe94c51b4ccbff7cb95840371f7d2d77b87ea9df527e8c9c',
  (3, 100): 'ba3f3d1e70fe73c666bede150ca7dbd0f332fc02959fe5178f8de8141b712b14',
  (5, 1):   '3a42dd193fad5edfacbbddb6642510a6b009e3f5e7787122eb05fc610580131d',
  (5, 10):  '9818562e125a03cd529a4eaf440e444387c7563855744302e35bdac9f85fae56',
  (5, 100): '43dd50d6f91d9d85085558e66e2ec0168b403ded47c6dad43cd2acfddca2f618'
}

def fixed_wt(n, w):
    """ Create a n-bit vector with Weight w. """
    x = 0
    flip = n > w // 2
    if flip:
        w = n - w
    for _ in range(w):
        i = random.randrange(n)
        while (x >> i) & 1:
            i = random.randrange(n)
        x |= 1 << i
    if flip:
        x = (~x) & ((1 << n) - 1)

    return x


#   Test Reed-Solomon with error weight w

def test_rs(iut, w):
    k   =   iut.k           #   message length
    t   =   2 * iut.delta   #   code length

    #   RS Encode a random message
    m0  = random.randbytes(k)
    #print('m0:', m0.hex())

    c0  = iut.rs_encode(m0)
    #print('c0:', c0.hex())

    #   random error wector with weight w
    e   = bytearray(t + k)
    i   = 0
    while i < w:
        j = random.randrange(t + k)
        if e[j] == 0:
            e[j] = random.randrange(255) + 1
            i += 1

    #print('e: ', e.hex())
    c1  =   c0.copy()
    for i in range(t + k):
        c1[i]   ^= e[i]
    #print('c1:', c1.hex())

    #   decode the vector with error
    m1  =   iut.rs_decode(c1)
    #print('m1:', m1.hex())

    ok  =   m0 == m1
    return ok

#   Test Reed-Muller

def test_rm(iut, w):

    #print(f'=== Reed-Muller test HQC-{iut.lev}  w = {w}')

    #   RS Encode a random message
    m0  = random.randbytes(iut.n1)
    #print('m0:', m0.hex())

    c0  = iut.rm_encode(m0)
    #print('c0:', hex(c0))

    e   = fixed_wt(iut.n1 * iut.n2, w)
    #print('e: ', hex(e))
    c1  = c0 ^ e
    #print('c1:', hex(c1))

    m1  = iut.rm_decode(c1)

    #   get the byte difference (as that matters for RS cascade.)
    d   = 0
    for i in range(iut.n1):
        if m0[i] != m1[i]:
            d += 1

    return d


#   HQC uses its own testing PRNG based on SHAKE256
def prng_init(entropy_input=b'', personalization_string=b''):
    prng = SHAKE256.new(entropy_input + personalization_string + b'\00')
    return prng

def nist_kat(lab='', x=b''):
    return f'{lab} = {x.hex().upper()}\n'

#   generate a string output equivalent to HQC's .rsp files
def hqc_print_rsp(iut, katnum=10):
    r       = f'# HQC-{iut.lev}\n\n'
    prng0   = prng_init(bytes(range(48)))
    print(f'[TEST] {iut.alg_id} KAT {katnum} ', end='', flush=True)
    for i in range(katnum):
        r       +=  f'count = {i}\n'
        seed0    =  prng0.read(48)
        r       +=  nist_kat('seed', seed0)
        prng     =  prng_init(seed0)
        (pk, sk) =  iut.keygen(prng)
        r       +=  nist_kat('pk', pk)
        r       +=  nist_kat('sk', sk)

        (ss, ct) =  iut.kem_encaps(prng, pk)
        r       +=  nist_kat('ct', ct)
        r       +=  nist_kat('ss', ss)

        s2       =  iut.kem_decaps(sk, ct)
        if ss == s2:
            print('.', end='', flush=True)
        else:
            print(f'\n[FAIL] {iut.alg_id} KAT decaps failure @ {i}')
            exit(0)
        r   +=      '\n'
    print(' [END]')
    return r


if (__name__ == "__main__"):

    #   test vectors
    katnum = 100
    for iut in HQC_ALL:
        kat =   hqc_print_rsp(iut, katnum)
        md  =   SHA256.new(kat.encode('ASCII')).hexdigest()
        if HQC_KAT_SHA256[(iut.lev, katnum)] != md:
            print(f'[FAIL] {iut.alg_id} KAT {katnum} {md}')
            exit(0)
        else:
            print(f'[PASS] {iut.alg_id} KAT {katnum} {md}')

    #   error correction
    for iut in HQC_ALL:
        p = 0.0
        r = 0.3
        while r <= 0.5:
            w = int(r * iut.n1n2)
            d = test_rm(iut, w)
            ok = d <= iut.delta
            print(f'RM {iut.alg_id}  p= {r:5.3f}  w= {w}  d= {d}  {ok}')
            if ok:
                p = r
                r += 0.01
            else:
                break

        #   0.4 is a bit arbitrary
        if p < 0.4:
            print(f'[FAIL] {iut.alg_id} RM p {p} too low.')
            exit(0)
        else:
            print(f'[PASS] {iut.alg_id} RM p {p:5f} ok.')

        w = 0
        i = 0
        while i <= iut.n1:
            ok = test_rs(iut, i)
            print (f'RS {iut.alg_id}  w= {i}  {ok}')
            if ok:
                w = i
                i += 1
            else:
                break
        if w < iut.delta:
            print(f'[FAIL] {iut.alg_id} RS dist {w} < {iut.delta}.')
            exit(0)
        else:
            print(f'[PASS] {iut.alg_id} RS dist {w} >= {iut.delta}.')


