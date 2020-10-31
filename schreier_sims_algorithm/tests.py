from permutation import *
from schreier_tree import *
from stabilizer_chain import *

# !!!АХТУНГ!!!
# Запускать этот тест только при серьезной необходимости
# Мой комплюктер пыхтел 10 минут (ответ благо правильный)
def hard_test():
    check = []
    for i in range(3, 55):
        check.append([1, 2, i])

    h = set()
    for e in check:
        h.add(Cycle(e).to_permutation())
    checker = StabilizerChain(h)

    ans = 1
    for x in range(2, 55):
        ans *= x
    if ans == checker.group_ord():
        print("Тест 99 уровня пройден. Легендарная победа")
    else:
        print("Тест 99 уровня провален")

def easy_peasy_lemon_squeezy():
    check = []
    for i in range(3, 6):
        check.append([1, 2, i])

    h = set()
    for e in check:
        h.add(Cycle(e).to_permutation())
    checker = StabilizerChain(h)

    print(checker.group_ord())
    print(checker.contains(Cycle([1, 2, 3]).to_permutation()))

def serious_test():
    perms = set()
    tmp1 = Cycle("(1 2 3)").to_permutation()
    tmp2 = Cycle("(1 2 4)").to_permutation()
    perms.add(tmp1)
    perms.add(tmp2)

    checker = StabilizerChain(perms)

    print(checker.group_ord())
    perm1 = Cycle("(1 2)").to_permutation()
    perm2 = Cycle("(3 4)").to_permutation()
    answer = checker.contains(perm1 * perm2)

    str = ""
    for e in answer[1]:
        str = str + e.__str__()
    print(str)