import random
import additional_funcs
import math


def miller_rabin(n):
    m = n - 1
    s = 0
    while m % 2 == 0:
        s += 1
        m //= 2

    for k in range(int(math.log(n, 2)) + 1):
        a = random.randint(1, n - 1)
        x = additional_funcs.fast_pow(a, m, n)

        for i in range(s):
            if ((x ** 2) % n) == 1 and x != 1 and x != n - 1:
                return False
            x = (x ** 2) % n

        if x != 1:
            return False

    return True


def gen_pseudo_primes():
    y = True
    while y:
        n = random.randint(2 ** 123, 2 ** 128)

        if miller_rabin(n):
            return n
