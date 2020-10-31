import random
import math
import additional_funcs

primes = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127]

def gen_num():
    n = 2
    unique = set()
    unique.add(n)
    while n < 2 ** 123:
        next_prime = random.randint(0, len(primes) - 1)
        n *= primes[next_prime]
        unique.add(primes[next_prime])

    if len(unique) < 5 or n < 2 ** 128:
        return gen_num()

    return sorted(list(unique)), n + 1


def gen_primes():
    while True:
        cur_primes, n = gen_num()
        a = luke_test(n, cur_primes)
        if a != False:
            return n, cur_primes, a


def luke_test(n, cur_primes):
    for a in range(2, int(math.log(n, 2)) + 1):
        if additional_funcs.fast_pow(a, n - 1, n) != 1:
            return False

        for p in cur_primes:
            if additional_funcs.fast_pow(a, (n - 1) // p, n) != 1:
                return a
