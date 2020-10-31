import rsa
import gen_primes

n, p, q, e, d = rsa.gen_keys()
print("Открытый ключ - (", "n = ", n, "e = ", e, ")")
print("Закрытый ключ - (p = ", p, "q = ", q, "d = ", d, ")")

print("Введите сообщение t = ")
t = int(input())

s = rsa.rsa_encrypt(n, e, t)
print("Зашифрованное сообщение s = ", s)

if rsa.rsa_decrypt(n, d, s) == t:
    print("Сообщение расшифрованно успешно!\n")