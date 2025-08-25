#   hqc.py
#   2025-03-31  Markku-Juhani O. Saarinen <mjos@iki.fi>

#   === Currently matches 2025-08-22 version of the spec.

from Crypto.Hash import SHAKE256, SHA3_256, SHA3_512

def dbg_dump(lab='', x=b''):
    l = len(x)
    print(f'{lab} = ({l})')
    for i in range(0, l, 32):
        print(f'{lab}[{i:06x}]: {x[i:i+32].hex()}')

def dbg_hex(lab='', x=b''):
    l = len(x)
    print(f'{lab}[{l}] =', x.hex().upper())

def dbg_chk(lab='', data=b''):
    def crc32byte(x, c):
        x ^= c << 24
        for _ in range(8):
            x <<= 1
            if x & 0x100000000:
                x ^= 0x104C11DB7
        return x
    x = 0
    for y in data:
        x = crc32byte(x, y)
    l = len(data)
    while l > 0:
        x = crc32byte(x, l & 0xFF)
        l >>= 8
    x = (~x) & 0xFFFFFFFF;
    print(f'{lab}: {x:08X} ({len(data)})')

class HQC:

    def __init__(self, lev, n, n1, n2, w, w_e, w_r, delta):
        self.alg_id     =   f'HQC-{lev}'            # identifier
        self.lev        =   lev                     # { 1, 3, 5 }
        self.sec_sz     =   (4 * self.lev) + 12     # { 16, 24, 32 }
        self.n          =   n                       # ring x^n + 1 size
        self.n_sz       =   (n + 7) // 8            # ring size in bytes
        self.n_rej      =   (2**24 // n) * n        # rejection bound
        self.n1         =   n1                      # Reed-Solomon code len
        self.n2         =   n2                      # Reed-Muller code len
        self.n1n2       =   n1 * n2                 # product code, n1*n2 < n
        self.n1n2_sz    =   (self.n1n2 + 7) // 8    # size in bytes
        self.w          =   w                       # x, y weight (keygen)
        self.w_e        =   w_e                     # e weight (encrypt)
        self.w_r        =   w_r                     # r1, r2 weight (encrypt)
        self.k          =   self.sec_sz             # encoded message
        self.delta      =   delta
        self.salt_sz    =   16                      # salt (ciphertext)
        self.seed_sz    =   32                      # general seed length

        #   external parameters
        self.pk_sz      =   self.seed_sz + self.n_sz
        self.sk_sz      =   self.pk_sz + self.seed_sz + self.k + self.seed_sz
        self.ct_sz      =   self.n_sz + self.n1n2_sz + self.salt_sz
        self.ss_sz      =   32

        self.alpha      =   0b10                    # generator in GF(256)
        self.rs         =   self.rs_genpoly(self.delta * 2)

    #   === Internal

    def xof_init(self, seed=b''):
        xof = SHAKE256.new(seed + b'\01')
        return xof

    def xof_get_bytes(self, xof, sz):
        """ this functionally matches with the *implementation* """
        b = xof.read(sz)
        if sz % 8 != 0:             #   need to waste remainder
            xof.read(8 - sz % 8)
        return b

    def sample_fixed_wt_mod(self, xof, wt):
        """ Random vector of given weight, using modular reduction, shuffle. """
        rand_u32 = self.xof_get_bytes(xof, 4 * wt)
        supp = []
        for i in range(wt):
            u32 = int.from_bytes(rand_u32[4*i : 4*(i+1)], byteorder='little')
            supp += [ ((u32 * (self.n - i)) >> 32) + i ]
            #print(f'[{i:2}] {supp[i]:5}')

        for i in range(wt - 1, -1, -1):
            for j in range(i + 1, wt):
                if supp[i] == supp[j]:
                    supp[i] = i
                    break
        v = 0
        for i in range(wt):
            v |= 1 << supp[i]
        return v

    def sample_fixed_wt_rej(self, ctx_pke_dk, wt):
        """ Random vector of given weight, using rejection sampling. """
        v = 0
        i = 0
        j = 3 * wt
        while i < wt:
            #   reference code uses 3 * wt chunks, we will have to do so too
            if j >= (3 * wt):
                b = self.xof_get_bytes(ctx_pke_dk, 3 * wt)
                j = 0
            #   reference code is big endian here
            sup = int.from_bytes(b[j:j+3], byteorder='big')
            j += 3

            if sup < self.n_rej:    #   n_rej: largest multiple of n < 2**24
                sup %= self.n
                bit = 1 << sup
                if (v & bit) == 0:
                    v |= bit
                    i += 1
        return  v

    def sample_vect(self, xof):
        """ Random vector """
        v = self.xof_get_bytes(xof, self.n_sz)
        v = int.from_bytes(v, byteorder='little')
        v &= (1 << self.n) - 1
        return v

    def vect_mul(self, a, b):
        """ Multiply binary polynomials: a.b  mod x^n+1 """
        r = 0
        for i in range(self.n):
            if (a >> i) & 1:
                r ^= b << i
        r = (r ^ (r >> self.n)) & ((1 << self.n) - 1)
        return r

    def gf_mul(self, a, b):
        """ GF(256): multiply a*b mod (x^8 + x^4 + x^3 + x^2 + 1). """
        r = a & (-(b & 1));
        for i in range(1, 8):
            a = (a << 1) ^ ((-(a >> 7)) & 0x11D);
            r ^= a & (-((b >> i) & 1));
        return r

    def gf_exp(self, a, e):
        """ GF(256): Compute a**e in the finite field. """
        r = 1
        if e & 1:
            r = a
        e >>= 1
        while e > 0:
            a = self.gf_mul(a, a)
            if e & 1:
                r = self.gf_mul(r, a)
            e >>= 1
        return r

    def gf_inv(self, x):
        """ GF(256): Multiplicative inverse 1/x. """
        return self.gf_exp(x, 254)

    def gf_gauss(self, m):
        """ Gaussian elimination: Create I on the left. """
        r   = len(m)
        c   = len(m[0])

        for i in range(r):
            j = i
            while j < r and m[j][i] == 0:
                j += 1
            if j >= r:
                continue
            if j > i:
                m[i],m[j] = m[j],m[i]
            x = self.gf_inv(m[i][i])
            for k in range(c):
                m[i][k] = self.gf_mul(x, m[i][k])
            for j in range(r):
                if i == j or m[j][i] == 0:
                    continue
                x = m[j][i]
                for k in range(c):
                    m[j][k] ^= self.gf_mul(x, m[i][k])
        return m

    #   === Coders

    def rs_genpoly(self, t):
        """ Reed-Solomon: Compute generator polynomial. """
        f   =   bytearray([1])      #   start with degree-0 polynomial 1
        a   =   1
        for i in range(t):
            a   = self.gf_mul(a, self.alpha)    #   a = alpha^(i+1)
            f.insert(0, 0)          #   multiply by x (shift left)
            for j in range(i+1):    #   multipliyng by (x + a)
                f[j] ^= self.gf_mul(a, f[j + 1])
        return f

    def rs_encode(self, msg):
        """ Reed-Solomon encoder. """
        y = bytearray(self.n1 - self.k)
        for i in range(self.k):
            x = msg[self.k - 1 - i] ^ y[-1]
            y = bytearray(1) + y[:-1]
            for j in range(self.n1 - self.k):
                y[j] ^= self.gf_mul(x, self.rs[j])
        return y + msg

    def rm_encode(self, msg):
        """ Reed-Muller RM(1,7). """
        HQC_RM_TAB  =   [   #   Reed-Muller matrix
            0xAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA,
            0xCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC,
            0xF0F0F0F0F0F0F0F0F0F0F0F0F0F0F0F0,
            0xFF00FF00FF00FF00FF00FF00FF00FF00,
            0xFFFF0000FFFF0000FFFF0000FFFF0000,
            0xFFFFFFFF00000000FFFFFFFF00000000,
            0xFFFFFFFFFFFFFFFF0000000000000000,
            0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
        ]
        r = 0
        p = 0
        rep = self.n2 // 128
        for x in msg:
            y = 0
            for i in range(8):
                if (x >> i) & 1:
                    y ^= HQC_RM_TAB[i]
            for _ in range(rep):
                r ^= y << p
                p += 128
        return r

    #   === Decoders

    def rm_decode(self, cw):
        """ Reed-Muller decode bytes from blocks. """
        m   = []
        #   for each byte
        for i in range(self.n1):

            #   sum the repeat codeblocks
            s = [0] * 128
            for j in range(0, self.n2, 128):
                for k in range(128):
                    s[k] += (cw >> ((i * self.n2) + j + k)) & 1

            #   Hadamard transform
            for j in range(7):
                s0, s1  =   [], []
                for k in range(0, 128, 2):
                    s0  +=  [ s[k] + s[k + 1] ]
                    s1  +=  [ s[k] - s[k + 1] ]
                s = s0 + s1

            #   subtract "constant term"
            s[0] -= 64 * (self.n2 // 128)

            #   find peaks ("Green machine")
            x, y, z = 0, -1, 0
            for j in range(128):
                if abs(s[j]) >  y:
                    x   = j
                    z   = s[x]
                    y   = abs(z)
            if  z > 0:
                x += 128
            m += [ x ]
        return bytearray(m)

    def rs_synd(self, cw, t = None):
        """ Reed-Solomon decode: compute syndrome. """
        ai  = 1
        if  t == None:
            t   = 2 * self.delta
        s   = bytearray(t)
        for i in range(t):
            ai  = self.gf_mul(ai, self.alpha)
            aij = ai
            x   = cw[0]                     #   alpha^0
            for j in range(1, self.n1):
                x   ^=  self.gf_mul(cw[j], aij)
                aij =   self.gf_mul(aij, ai)
            s[i] = x
        return s

    def rs_elp(self, syn):
        """ Reed-Solomon decode: Error locator polynomial. """
        d_sig   =   0   #   degrees
        d_sig_p =   0
        d_sig_t =   0
        x_sig_p =   bytearray(b'\0\1') + bytearray(self.delta)
        sig     =   bytearray(b'\1') + bytearray(self.delta)

        pp      =   -1
        dp      =   1

        for mu in range(2 * self.delta):
            d = syn[mu]
            for i in range(min(mu, self.delta)):
                d ^= self.gf_mul(sig[i + 1], syn[mu - i - 1])
            sig_t =  sig[:self.delta].copy()
            d_sig_t = d_sig
            dd  = self.gf_mul(d, self.gf_inv(dp))
            for i in range(min(mu + 1, self.delta)):
                sig[i + 1] ^= self.gf_mul(dd, x_sig_p[i + 1])
            d_x = mu - pp
            d_x_sig_p = d_x + d_sig_p
            if d != 0 and d_x_sig_p > d_sig:
                d_sig = d_x_sig_p
                pp  = mu
                dp  = d
                x_sig_p = b'\0' + sig_t
                d_sig_p = d_sig_t
            else:
                x_sig_p = b'\0' + x_sig_p[:-1]

        return sig

    def rs_roots(self, sig):
        """ Reed-Solomon decode: Get error locations. """
        l   = []
        deg = len(sig) - 1
        while deg > 0 and sig[deg] == 0:
            deg -= 1
        if deg == 0:
            return l
        x = 1
        g = self.gf_inv(2)
        for a in range(self.n1):
            y = sig[deg]
            for i in range(deg - 1, -1, -1):
                y = self.gf_mul(y, x)
                y ^= sig[i]
            if y == 0:
                l += [a]
            x = self.gf_mul(x, g)
        return l

    def rs_decode(self, cw):
        """ Reed-Solomon decode (slow). """
        syn = self.rs_synd(cw)      #   syndrome
        sig = self.rs_elp(syn)      #   error location polynomial
        ep  = self.rs_roots(sig)    #   error positions
        ne  = len(ep)               #   number of errors
        if ne == 0:
            return cw[-self.k:]     #   return early if none

        sm  = []
        for i in range(ne):         #   create simultaneous equations
            v = bytearray(len(cw))
            v[ep[i]] = 1
            u = bytearray(ne)
            u[i] = 1
            sm += [ self.rs_synd(v, ne) + u ]

        self.gf_gauss(sm)           #   Gaussian elimination

        ev  =   bytearray(ne)       #   substitute
        for i in range(ne):
            x = syn[i]
            for k in range(ne):
                ev[k] ^= self.gf_mul(x, sm[i][ne + k])

        for i in range(ne):         #   correct the errors
            cw[ep[i]] ^= ev[i]

        m = cw[-self.k:]            #   return the message part
        return m

    def pke_keygen(self, seed_pke):
        """ HQC-PKE.Keygen(seedPKE) """

        #   Compute dkPKE and ekPKE seeds
        #   1. (seedPKE.dk, seedPKE.ek) <- I(seedPKE)
        i_res       = SHA3_512.new(seed_pke + b'\02').digest()
        seed_pke_dk = i_res[0:32]
        seed_pke_ek = i_res[32:64]

        #   Compute decryption key dkPKE
        #   2.  ctxPKE.dk <- XOF.Init(seedPKE.dk)
        ctx_pke_dk  = self.xof_init(seed_pke_dk)

        #   3.  (ctxPKE.dk, y) <- SampleFixedWeightVect_$(ctxPKE.dk, Rw)
        y   = self.sample_fixed_wt_rej(ctx_pke_dk, self.w)

        #   4.  (ctxPKE.dk, x) <- SampleFixedWeightVect_$(ctxPKE.dk, Rw)
        x   = self.sample_fixed_wt_rej(ctx_pke_dk, self.w)

        #   5.  dkPKE <- seedPKE.dk
        dk_pke  =   seed_pke_dk

        #   Compute encryption key ekPKE
        #   6.  ctxPKE.ek <- XOF.Init(seedPKE.ek)
        ctx_pke_ek  = self.xof_init(seed_pke_ek)

        #   7.  (ctxPKE.ek, h) <- SampleVect(ctxPKE.ek, R)
        h   =   self.sample_vect(ctx_pke_ek)

        #   8.  s <- x + h * y
        s   =   self.vect_mul(y, h) ^ x

        #   9.  ekPKE <- (seedPKE.ek, s)
        ek_pke  =   seed_pke_ek + s.to_bytes(self.n_sz, byteorder='little')

        #   10. return (ekPKE, dkPKE)
        return (ek_pke, dk_pke)

    def keygen(self, xkg):
        """ HQC-KEM.Keygen() """

        #   Sample seedKEM
        #   1. seedKEM <-$ B|seed|
        seed_kem    = xkg.read(self.seed_sz)

        #   Compute seedPKE and randomness sigma
        #   2.  ctxKEM <- XOF.Init(seedKEM)
        ctx_kem     = self.xof_init(seed_kem)

        #   3.  (ctxKEM, seedPKE) <- XOF.GetBytes(ctxKEM, |seed|)
        seed_pke    = ctx_kem.read(self.seed_sz)

        #   4.  (ctxKEM, sigma) <- XOF.GetBytes(ctxKEM, |k|)
        sigma       = ctx_kem.read(self.sec_sz)

        #   Compute HQC-PKE keypair
        #   5.  (ekPKE, dkPKE) <- HQC-PKE.Keygen(seedPKE)
        (ek_pke, dk_pke) = self.pke_keygen(seed_pke)

        #   Compute HQC-KEM keypair
        #   6.  ekKEM <- ekPKE
        ek_kem  = ek_pke

        #   7.  dkKEM <- (ekKEM, dkPKE, sigma, seedKEM)
        dk_kem  = ek_kem + dk_pke + sigma + seed_kem

        #   8.  return (ekKEM, dkKEM)
        return  (ek_kem, dk_kem)


    def pke_encrypt(self, ek_pke, m, theta):
        """ HQC-PKE.Encrypt(ekPKE, m, theta) """

        #   Parse encryption key ekPKE
        #   1.  seedPKE.ek <- ekPKE [0 : |seed|]
        seed_pke_ek =   ek_pke[0 : self.seed_sz]

        #   2.  ctxPKE.ek <- XOF.Init(seedPKE.ek)
        ctx_pke_ek  = self.xof_init(seed_pke_ek)

        #   3.  (ctxPKE.ek, h) <- SampleVect(ctxPKE.ek, R)
        h   =   self.sample_vect(ctx_pke_ek)

        #   4.  s <- ekPKE [|seed| : |seed| + |s|]
        s   =   int.from_bytes(
                    ek_pke[self.seed_sz : self.seed_sz + self.n_sz],
                        byteorder='little')

        #   Compute ciphertext cPKE
        #   5.  ctx_th <- XOF.Init(theta)
        ctx_th  =   self.xof_init(theta)

        #   6.  (ctx_th, r2) <- SampleFixedWeightVect(ctx_th, Rw_r)
        r2  = self.sample_fixed_wt_mod(ctx_th, self.w_r)

        #   7.  (ctx_th, e) <- SampleFixedWeightVect(ctx_th, Rw_e)
        e   = self.sample_fixed_wt_mod(ctx_th, self.w_e)

        #   8.  (ctx_th, r1) <- SampleFixedWeightVect(ctx_th, Rw_r)
        r1  = self.sample_fixed_wt_mod(ctx_th, self.w_r)

        #   9.  u <- r1 + h * r2
        u   = self.vect_mul(r2, h) ^ r1
        u   = u.to_bytes(self.n_sz, byteorder='little')

        #   10. v <- C.Encode(m) + Truncate(s * r2 + e, l)
        cm  = self.rm_encode(self.rs_encode(m))     #   RM(DS(m))
        v   = cm ^ self.vect_mul(r2, s) ^ e
        v   &= (1 << self.n1n2) - 1
        v   = v.to_bytes(self.n1n2_sz, byteorder='little')

        #   11. cPKE <- (u, v)
        c_pke   =   u + v

        #   12. return cPKE
        return  c_pke

    def kem_encaps(self, prng, ek_kem):
        """ HQC-KEM.Encaps(ekKEM) """

        #   Sample message m and salt
        #   1.  m <-$ B|k|
        m       = prng.read(self.k)

        #   2.  salt <-$ B|salt|
        salt    = prng.read(self.salt_sz)

        #   Compute shared key K and ciphertext cKEM
        #   3.  (K, theta) <- G(H(ekKEM) || m || salt)
        tmp_h   = SHA3_256.new(ek_kem + b'\01').digest()
        tmp_g   = SHA3_512.new(tmp_h + m + salt + b'\00').digest()
        kk      = tmp_g[0:32]
        theta   = tmp_g[32:64]

        #   4.  cPKE <- HQC-PKE.Encrypt(ekKEM, m, theta)
        c_pke   = self.pke_encrypt(ek_kem, m, theta)

        #   5.  cKEM <- (cPKE, salt)
        c_kem   = c_pke + salt

        #   6.  return (K, cKEM)
        return  (kk, c_kem)


    def pke_decrypt(self, dk_pke, c_pke):
        """ HQC-PKE.Decrypt(dkPKE, cPKE) """

        #   Parse decryption key dkPKE
        #   1.  seedPKE.dk <- dkPKE[0 : |seed|]
        seed_pke_dk = dk_pke[0 : self.seed_sz]

        #   2.  ctxPKE.dk <- XOF.Init(seedPKE.dk)
        ctx_pke_dk  = self.xof_init(seed_pke_dk)

        #   3.  (ctxPKE.dk, y) <- SampleFixedWeightVect$(ctxPKE.dk, Rw)
        y   = self.sample_fixed_wt_rej(ctx_pke_dk, self.w)

        #   Parse ciphertext cPKE
        #   4.  u <- cPKE[0 : |u|]
        u   = int.from_bytes(c_pke[0 : self.n_sz], byteorder='little')

        #   5.  v <- cPKE[|u| : |u| + |v|]
        v   = int.from_bytes(c_pke[self.n_sz : self.n_sz + self.n1n2_sz],
                                byteorder='little')

        #   Compute plaintext m
        #   6.  m <- C.Decode(v - Truncate(u * y, l))
        cm  = v ^ self.vect_mul(y, u)
        m   = self.rm_decode(cm)
        m   = self.rs_decode(m)

        #   7.  return m
        return  m

    #   === KEM functionality

    def kem_decaps(self, dk_kem, c_kem):
        """ HQC-KEM.Decaps(dkKEM , cKEM) """

        #   Parse decapsulation key dkKEM
        #   1.  ekKEM <- dkKEM[0 : |ekKEM|]
        ek_kem  =   dk_kem[0 : self.pk_sz]

        #   2.  dkPKE <- dkKEM[|ekKEM| : |ekKEM| + |dkPKE|]
        dk_pke  =   dk_kem[self.pk_sz : self.pk_sz + self.seed_sz]

        #   3.  sigma <- dkKEM[|ekKEM| + |dkPKE| : |ekKEM| + |dkPKE| + |sigma|]
        sigma   =   dk_kem[self.pk_sz + self.seed_sz :
                            self.pk_sz + self.seed_sz + self.sec_sz]

        #   Parse ciphertext cKEM
        #   4.  cPKE <- cKEM[0 : |cPKE|]
        c_pke   =   c_kem[0 : self.n_sz + self.n1n2_sz]

        #   5.  salt <- cKEM[|cPKE| : |cPKE| + |salt|]
        salt    =   c_kem[self.n_sz + self.n1n2_sz : self.ct_sz]

        #   Compute message m'
        #   6.  m' <- HQC-PKE.Decrypt(dkPKE, cPKE)
        m_p     =   self.pke_decrypt(dk_pke, c_pke)

        #   Compute shared key K' and ciphertext c'KEM
        #   7.  (K', theta') <- G(H(ekKEM) || m' || salt)
        tmp_h   = SHA3_256.new(ek_kem + b'\01').digest()
        tmp_g   = SHA3_512.new(tmp_h + m_p + salt + b'\00').digest()
        kk_p    = tmp_g[0:32]
        theta_p = tmp_g[32:64]

        #   8.  c'PKE <- HQC-PKE.Encrypt(ekKEM, m', theta')
        c_pke_p = self.pke_encrypt(ek_kem, m_p, theta_p)

        #   9.  c'KEM <- (c'PKE, salt)
        c_kem_p = c_pke_p + salt

        #   Compute rejection key K~
        #   10. K~ <- J(H(ekKEM) || sigma || cKEM)
        k_rej   = SHA3_256.new(tmp_h + sigma + c_kem + b'\03').digest()

        #   11. if m' = fail or if c'KEM != cKEM
        #   12.     K' <- K~
        #   13. endif
        if m_p == None or c_kem_p != c_kem:
            kk_p =  k_rej

        #   14. return K'
        return  kk_p

#   === Parameter sets

HQC_1   = HQC(  lev = 1,    n = 17669,  n1  = 46,   n2  = 384,
                w   = 66,   w_e = 75,   w_r = 75,   delta = 15  )
HQC_3   = HQC(  lev = 3,    n = 35851,  n1  = 56,   n2  = 640,
                w   = 100,  w_e = 114,  w_r = 114,  delta = 16  )
HQC_5   = HQC(  lev = 5,    n = 57637,  n1  = 90,   n2  = 640,
                w   = 131,  w_e = 149,  w_r = 149,  delta = 29  )

HQC_ALL = [ HQC_1, HQC_3, HQC_5 ]

#   === KAT Generator fot testing purposes

if (__name__ == "__main__"):

    #   HQC uses its own testing PRNG based on SHAKE256
    def prng_init(entropy_input=b'', personalization_string=b''):
        prng = SHAKE256.new(entropy_input + personalization_string + b'\00')
        return prng

    def nist_kat(lab='', x=b''):
        print(f'{lab} = {x.hex().upper()}')

    #   generate std output equivalent to HQC's .rsp files
    def hqc_print_rsp(iut):
        print(f'# HQC-{iut.lev}\n')
        prng0   = prng_init(bytes(range(48)))
        for i in range(100):
            print(f'count = {i}')
            seed0   = prng0.read(48)
            nist_kat('seed', seed0)
            prng    = prng_init(seed0)

            (pk, sk) = iut.keygen(prng)
            nist_kat('pk', pk)
            nist_kat('sk', sk)

            (ss, ct) = iut.kem_encaps(prng, pk)
            nist_kat('ct', ct)
            nist_kat('ss', ss)

            s2 = iut.kem_decaps(sk, ct)
            if s2 != ss:
                print('ERROR -- decryption failure')
            print()

            c2  = bytearray(ct)
            c2[0] ^= 1

            s2 = iut.kem_decaps(sk, c2)
            dbg_hex('s2', s2)

    #   our test case
    hqc_print_rsp(HQC_1)


