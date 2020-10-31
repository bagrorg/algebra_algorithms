def fast_pow(m, z, n):
    if z == 0:
        return 1
    c = fast_pow(m, z // 2, n)
    if z % 2 == 0:
        return (c ** 2) % n
    else:
        return (m * (c ** 2)) % n


def gcd_extended(a, b):
    x0, x1, y0, y1 = 0, 1, 1, 0
    while a != 0:
        q, b, a = b // a, a, b % a
        y0, y1 = y1, y0 - q * y1
        x0, x1 = x1, x0 - q * x1
    return x0, y0


def get_opposite(e, phi):
    x, y = gcd_extended(e, phi)
    if x < 0:
        x += phi
    return x