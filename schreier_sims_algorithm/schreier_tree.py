from permutation import *


class SchreierTree:
    def __init__(self, w, S, n):
        # Создаем корень
        self.root = SchrierNode(w)

        # Орбита будет иметь тип Словарь,
        # так как так легче с ней обращаться.
        # в orbit[s] будет лежать такой g, что gw = s
        orbit = dict()
        
        # Для w это g = id 
        orbit[w] = Permutation(n + 1)
        
        # Обходом в глубину создадим дерево Шрайера
        self.dfs_build(w, S, orbit, self.root)

        self.orbit = orbit

    def dfs_build(self, w, S, orbit, node):
        # По всем элементам из S
        for g in S:
            # Если этого элемента еще не было в орбите
            if g[w] not in orbit:
                # Для нового элемента перестановка будет
                # равняться g * h, где h: hw = sigma и g * sigma = g[w]
                orbit[g[w]] = g * orbit[w]

                # Запускаем рекурсию от новой вершины, в которой лежит новый
                # элемент и g^-1, при применении которой мы получим исходный
                # на данной итерации элемент, и от нового элемента
                new_node = self.dfs_build(g[w], S, orbit, SchrierNode(g[w], g ** (-1)))

                # Добавим нового соседа текущей вершине
                node.add_neigh(new_node)
        return node

    def dfs_return_edge_elements(self, root, res):
        if root.elem != self.root.elem:
            res.add(root.parent_perm)
        for e in root.neigh_set:
            res = self.dfs_return_edge_elements(e, res)
        return res

    def get_orbit(self):
        return self.orbit

    def __len__(self):
        return len(self.orbit)


class SchrierNode:
    def __init__(self, elem, perm=0):
        self.elem = elem

        self.parent_perm = perm
        self.neigh_set = []

    def add_neigh(self, other):
        self.neigh_set.append(other)
