#   unit_tests.py
#   2025-04-03  Markku-Juhani O. Saarinen <mjos@iki.fi>

#   === Unit test for error correction

from Crypto.Hash import SHAKE256,SHA256
import random
from hqc import HQC_ALL

#   sha-256 hashes of KAT rsp files

HQC_KAT_SHA256 = {
    ( 128,  1 ):
        '97997017316d249be02665a51530c904041c584f8a8cbfe126e71a0ca5095255',
    ( 128,  10 ):
        'e2f87a16ece9153c3d3b8d4fcce735c9ff4a920e73dc26530cad3d96b7e2c7f7',
    ( 128,  100 ):
        '2be3afb1efb98ce58d719e19824f1a1fcb53fedb6ca8e16bc34afa53ac0528c0',
    ( 192,  1 ):
        '952596ee4239e4ef3d3316d738ab457236bc72de956e0262511bab40c4d941c6',
    ( 192,  10 ):
        '0d2717d6d4fb1a2dd762b5fe5b50f00deab6329fde27d39c7ec9f371bb1c6817',
    ( 192,  100 ):
        'a77dd028b01a6554eebe5530fbccf4891e788c50a7efae9156e7ebc01bbea20a',
    ( 256,  1 ):
        '0bd4ae5862f42992d2a20a82aa46cd06ee14addbf71048b046070635e678224e',
    ( 256,  10 ):
        '317bcab981dd37739800ccb0ed08f42fd32eed34e7065524bf97a3c5eb4a1664',
    ( 256,  100 ):
        'a63d941b042b3992dc9d19f9b901548bfcd7e5b9cb92bb87724ccfbbabee8917'
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

    #print(f'=== Reed-Muller test HQC-{iut.sec}  w = {w}')

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
    prng = SHAKE256.new(entropy_input + personalization_string + b'\01')
    return prng

def nist_kat(lab='', x=b''):
    return f'{lab} = {x.hex().upper()}\n'

#   generate a string output equivalent to HQC's .rsp files
def hqc_print_rsp(iut, katnum=10):
    r       =   f'# HQC-{iut.sec}\n\n'
    prng0   =   prng_init(bytes(range(48)))
    print(f'[TEST] {iut.alg_id} KAT {katnum} ', end='', flush=True)
    for i in range(katnum):
        r   +=      f'count = {i}\n'
        seed0 =     prng0.read(48)
        r   +=      nist_kat('seed', seed0)
        prng =      prng_init(seed0)
        (pk, sk) = iut.keygen(prng)
        r   +=      nist_kat('pk', pk)
        r   +=      nist_kat('sk', sk)

        (ct, ss) = iut.kem_enc(prng, pk)
        r   +=      nist_kat('ct', ct)
        r   +=      nist_kat('ss', ss)

        s2 = iut.kem_dec(ct, sk)
        if ss == s2:
            print('.', end='', flush=True)
        else:
            print(f'[FAIL] {iut.alg_id} KAT decaps failure count = {i}')
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
        if HQC_KAT_SHA256[(iut.sec, katnum)] != md:
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


