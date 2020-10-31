from permutation import *
from schreier_tree import *


class StabilizerChain:
    def __init__(self, S):
        # Чтобы понимать, чему в Sn равно n
        self.univ_len = max([len(x) for x in S])

        # Приводим все перестановки к одной длине
        for g in S:
            g.expand(self.univ_len)

        # Создаем цепочку
        self.chain = []

        # Инициализируем G0
        self.chain.append(Stabilizer(S, 0, self.univ_len))

        # Построение цепочки
        for i in range(1, self.univ_len + 1):
            # По Лемме Шрайера строим Si на основе порождающего множества для Gi-1
            # и орбите элемента, лежащего в прошлом Gi-1, т.е. i.
            new_generating_set = self.make_gen(self.chain[i - 1].generating_set, self.chain[i - 1].tree.orbit)
            # Пронормализуем, чтобы не было экспоненциального роста
            final_generating_set = list(self.normalize(new_generating_set))

            # Добавляем стабилизатор в цепочку
            self.chain.append(Stabilizer(final_generating_set, i, self.univ_len))

            # Если получившийся набор = {e}, то заканчиваем
            if len(final_generating_set) == 0:
                self.base = [p for p in range(1, i + 1)]
                break

        # Создаем сильное порождающее множество
        self.strong_generating_set = set()
        self.strong_generating_set = self.stong_gen_make()

    # Создаем порожд. множество с помощью Леммы Шрайера
    def make_gen(self, S, orbit):
        G = set()
        for s in S:
            for u in orbit:
                tmp1 = s * orbit[u]
                tmp2 = orbit[s[u]] ** (-1)
                G.add(tmp2 * tmp1)
        return G

    def base(self):
        return self.base

    # Нормализуем множество, чтобы не было сильного роста
    def normalize(self, S):
        S_new = set()
        base = [{} for i in range(self.univ_len)]
        for s in S:
            for x in range(1, self.univ_len + 1):
                if s[x] != x:
                    if s[x] in base[x]:
                        s = (s ** (-1)) * base[x][s[x]]
                    else:
                        base[x][s[x]] = s
                        S_new.add(s)
                        break
        return S_new

    # Создаем сильное порождающее множество
    def stong_gen_make(self):
        res = set()
        for stab in self.chain:
            tmp = set()
            tmp = stab.tree.dfs_return_edge_elements(stab.tree.root, tmp)
            res = res | tmp
        return res

    def contains(self, sigma):
        decompose = []
        sigma.expand(self.univ_len)
        s = Permutation(self.univ_len)
        s.perm = sigma.perm[::]
        s.length = sigma.length

        for stab in self.chain:
            w = s[stab.index + 1]
            if w not in stab.tree.orbit:
                return False, []

            g = stab.tree.orbit[w] ** (-1)
            s = g * s

            if g != Permutation(self.univ_len):
                decompose.append(g)

        return True, decompose

    def group_ord(self):
        order = 1
        for step in self.chain:
            order *= len(step.tree)
        return order


class Stabilizer:
    def __init__(self, generating_set, index, length):
        self.generating_set = set(list(generating_set)[::])
        self.index = index

        if len(generating_set) == 0:
            self.generating_set.add(Permutation(length + 1))

        self.tree = SchreierTree(index + 1, generating_set, length)
