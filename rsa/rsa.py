import additional_funcs
import gen_pseudo_primes


def gen_keys():
    q = gen_pseudo_primes.gen_pseudo_primes()
    p = gen_pseudo_primes.gen_pseudo_primes()
    n = p*q
    phi = (p-1)*(q-1)
    e = gen_pseudo_primes.gen_pseudo_primes()

    d = additional_funcs.get_opposite(e, phi)

    return n, p, q, e, d



def rsa_encrypt(n, e, t):
    return additional_funcs.fast_pow(t, e, n)



def rsa_decrypt(n, d, s):
    return additional_funcs.fast_pow(s, d, n)
