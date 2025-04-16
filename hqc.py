#   hqc.py
#   2025-03-31  Markku-Juhani O. Saarinen <mjos@iki.fi>

#   === Currently matches 2025-02-19 version of the spec.

from Crypto.Hash import SHAKE256

class HQC:

    def __init__(self, sec, n, n1, n2, w, w_e, w_r, delta):
        self.alg_id     =   f'hqc-{sec}'
        self.sec        =   sec
        self.sec_sz     =   self.sec // 8
        self.n          =   n
        self.n_sz       =   (self.n + 7) // 8
        self.n1         =   n1
        self.n2         =   n2
        self.n1n2       =   n1 * n2
        self.n1n2_sz    =   (self.n1n2 + 7) // 8
        self.w          =   w
        self.w_e        =   w_e
        self.w_r        =   w_r
        self.k          =   self.sec_sz
        self.delta      =   delta
        self.salt_sz    =   16
        self.seed_sz    =   40

        self.pk_sz      =   self.seed_sz + self.n_sz
        self.sk_sz      =   self.seed_sz + self.sec_sz + self.pk_sz
        self.ct_sz      =   self.n_sz + self.n1n2_sz + self.salt_sz
        self.ss_sz      =   64

        self.alpha      =   0b10        # generator in GF(256)
        self.rs         =   self.rs_genpoly(self.delta * 2)

    #   === Internal

    def xof_init(self, seed=b''):
        xof = SHAKE256.new(seed + b'\02')
        return xof

    def xof_bytes(self, xof, l):
        x = xof.read((l + 7) & ~7)
        return x[:l]

    def vect_fixed_wt(self, xof, wt):
        """ Random vector with wt 1 bits. """
        rand_u32 = self.xof_bytes(xof, 4 * wt)
        supp = []
        for i in range(wt):
            u32 = int.from_bytes(rand_u32[4*i:4*(i+1)], byteorder='little')
            supp += [ u32 % (self.n - i) + i ]
        for i in range(wt - 1, -1, -1):
            for j in range(i + 1, wt):
                if supp[i] == supp[j]:
                    supp[i] = i
                    break
        v = 0
        for i in range(wt):
            v |= 1 << supp[i]
        return v

    def vect_random(self, xof):
        """ Random vector """
        v = xof.read(8 * ((self.n + 63) // 64))
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

        pp      =   -1      # 2*rho
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

    #   === CCA "PKE" functions (internal)

    def keygen(self, xkg):
        """ Keypair generation for both PKE and KEM """

        #   get random quantities
        sk_seed = xkg.read(self.seed_sz)
        sigma   = xkg.read(self.sec_sz)
        pk_seed = xkg.read(self.seed_sz)

        #   compute secret key
        xsk = self.xof_init(sk_seed)
        y   = self.vect_fixed_wt(xsk, self.w)
        x   = self.vect_fixed_wt(xsk, self.w)

        #   compute public key
        xpk = self.xof_init(pk_seed)
        h   = self.vect_random(xpk);
        s   = self.vect_mul(y, h) ^ x

        #   serialize keypair
        pk  = pk_seed + s.to_bytes(self.n_sz, byteorder='little')
        sk  = sk_seed + sigma + pk
        return  (pk, sk)

    def pke_encrypt(self, m, theta, pk):
        """ PKE Encrypt. """

        #   nested code  RM(DS(m))
        cm  = self.rm_encode(self.rs_encode(m))

        #   decode (h, s) from public key
        xpk = self.xof_init(pk[:self.seed_sz])
        h   = self.vect_random(xpk);
        s   = int.from_bytes(pk[self.seed_sz:], byteorder='little')

        #   generate (r1, r2, e)
        xen = self.xof_init(theta)
        r2  = self.vect_fixed_wt(xen, self.w_r)
        e   = self.vect_fixed_wt(xen, self.w_e)
        r1  = self.vect_fixed_wt(xen, self.w_r)

        #   u = r2 * h + r1
        u   = self.vect_mul(r2, h) ^ r1

        #   v = cm + s * r2 + e
        v   = cm ^ self.vect_mul(r2, s) ^ e
        v   &= (1 << self.n1n2) - 1

        #   encode into bytes
        u   = u.to_bytes(self.n_sz, byteorder='little')
        v   = v.to_bytes(self.n1n2_sz, byteorder='little')
        return  (u, v)

    def pke_decrypt(self, u, v, sk_seed):
        """ PKE Decrypt. """
        #   ciphertext to integers
        u   = int.from_bytes(u, byteorder='little')
        v   = int.from_bytes(v, byteorder='little')

        #   decode secret key
        xsk = self.xof_init(sk_seed[:self.seed_sz])
        y   = self.vect_fixed_wt(xsk, self.w)

        #   Compute v - u.y
        cm  = v ^ self.vect_mul(y, u)

        #   Decode
        m       = self.rm_decode(cm)
        m       = self.rs_decode(m)
        return  m

    #   === KEM functionality

    def kem_enc(self, prng, pk):
        m       = prng.read(self.k)
        #print('m=', m.hex())
        salt    = prng.read(self.salt_sz)

        #   theta
        tmp     = m + pk[:self.salt_sz * 2] + salt + b'\03'
        theta   = SHAKE256.new(tmp).read(self.seed_sz)
        (u, v)  = self.pke_encrypt(m, theta, pk)

        #   shared secret
        tmp     = m + u + v + b'\04'
        ss      = SHAKE256.new(tmp).read(64)

        #   ciphertext
        ct      = u + v + salt

        return  (ct, ss)

    def kem_dec(self, ct, sk):

        #   unpack ct
        u       = ct[:self.n_sz]
        v       = ct[self.n_sz:self.n_sz + self.n1n2_sz]
        salt    = ct[self.n_sz + self.n1n2_sz:self.ct_sz]

        #   unpack sk
        sk_seed = sk[:self.seed_sz]
        sigma   = sk[self.seed_sz : self.seed_sz + self.sec_sz]
        pk      = sk[self.seed_sz + self.sec_sz : self.sk_sz]
        #dbg_chk('pk from sk', pk)

        m       = self.pke_decrypt(u, v, sk_seed)

        #   Fujisaki-Okamoto re-encryption
        #   theta
        tmp     = m + pk[:self.salt_sz * 2] + salt + b'\03'
        theta   = SHAKE256.new(tmp).read(self.seed_sz)
        (u2,v2) = self.pke_encrypt(m, theta, pk)

        #print("m' =", m.hex())
        #print("sigma' =", sigma.hex())
        if  u != u2 or v != v2:
            m   = sigma

        #   shared secret
        tmp     = m + u + v + b'\04'
        ss      = SHAKE256.new(tmp).read(64)
        return  ss

#   === Parameter sets

HQC_128 = HQC(  sec = 128,  n = 17669,  n1  = 46,   n2  = 384,
                w   = 66,   w_e = 75,   w_r = 75,   delta = 15 )
HQC_192 = HQC(  sec = 192,  n = 35851,  n1  = 56,   n2  = 640,
                w   = 100,  w_e = 114,  w_r = 114,  delta = 16 )
HQC_256 = HQC(  sec = 256,  n = 57637,  n1  = 90,   n2  = 640,
                w   = 131,  w_e = 149,  w_r = 149,  delta = 29 )

HQC_ALL = [ HQC_128, HQC_192, HQC_256 ]

#   === KAT Generator fot testing purposes

if (__name__ == "__main__"):

    #   HQC uses its own testing PRNG based on SHAKE256
    def prng_init(entropy_input=b'', personalization_string=b''):
        prng = SHAKE256.new(entropy_input + personalization_string + b'\01')
        return prng

    def nist_kat(lab='', x=b''):
            print(f'{lab} = {x.hex().upper()}')

    #   generate std output equivalent to HQC's .rsp files
    def hqc_print_rsp(iut):
        print(f'# HQC-{iut.sec}\n')
        prng0   = prng_init(bytes(range(48)))
        for i in range(100):
            print(f'count = {i}')
            seed0   = prng0.read(48)
            nist_kat('seed', seed0)
            prng    = prng_init(seed0)

            (pk, sk) = iut.keygen(prng)
            nist_kat('pk', pk)
            nist_kat('sk', sk)

            (ct, ss) = iut.kem_enc(prng, pk)
            nist_kat('ct', ct)
            nist_kat('ss', ss)

            s2 = iut.kem_dec(ct, sk)
            if s2 != ss:
                print('ERROR')
            print()

    #   our test case
    hqc_print_rsp(HQC_128)

