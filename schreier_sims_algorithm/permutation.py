class Permutation:
    def __init__(self, n):
        # Не учитываем None
        self.length = n - 1
        self.perm = []
        self.perm.append(None)
        for x in range(1, n):
            # Создаем id
            self.perm.append(x)

    def __mul__(self, other):
        # Максимальная длина, чтобы не ошибиться с универсумом
        length = max(other.length, self.length)
        self.expand(length)
        other.expand(length)
        # +1 так как надо учесть None
        tmp = Permutation(length + 1)

        # Записываем в tmp нашу self перестановку
        for i in range(1, length + 1):
            tmp.perm[i] = other.perm[i]

        # Применяем перестановку other к tmp
        for x in range(1, other.length + 1):
            tmp.perm[x] = self.perm[tmp.perm[x]]
        return tmp

    def ord(self):
        # Создадим id, чтобы сравнить с ней. Создаем, чтобы не вызывать конструктор каждый раз в цикле
        id = Permutation(self.length + 1)

        # Создаем квадрат данной перестановки. Квдрат, так как легче обойти то, что при = self передается ссылка
        tmp = self * self
        power = 2
        while id != tmp:
            tmp = tmp * self
            power += 1
        return power

    def __pow__(self, power, modulo=None):
        # Создаем копию.
        # К сожалению =self в питоне передаст ссылку и мы будем менять исходный объект
        tmp = Permutation(self.length + 1)
        tmp.perm = self.perm[::]

        # Учтем отрицательную степень.
        # При power == -1 выдаст обратную перестановку.
        if power < 0:
            for i in range(1, self.length + 1):
                tmp.perm[self.perm[i]] = i

        for i in range(1, abs(power)):
            tmp = tmp * tmp
        return tmp

    def __eq__(self, other):
        length = max(self.length, other.length)
        self.expand(length)
        other.expand(length)
        for x in range(1, self.length + 1):
            if self.perm[x] != other.perm[x]:
                return False
        return True

    def __getitem__(self, item):
        return self.perm[item]

    def __len__(self):
        return self.length

    def to_cycles(self):
        # Массив циклов
        res = []

        # Чтобы не вызывать in каждый раз,
        # для скорости будем поддерживать used
        used = [False] * (self.length + 1)
        tmp = Cycle()

        for i in range(1, self.length + 1):
            # Если данный x \in X был рассмотрен, то пропускаем
            if used[i]:
                continue

            used[i] = True
            tmp.data.append(i)
            j = self[i]

            # Строим независимый цикл
            while j != tmp.data[1]:
                used[j] = True
                tmp.data.append(j)
                j = self[j]

            res.append(tmp)
            tmp = Cycle()

        return res

    # Для set() и dict()
    def __hash__(self):
        return hash(tuple(self.perm))

    def expand(self, length):
        while self.length < length:
            self.length += 1
            self.perm.append(self.length)

    def __str__(self):
        tmp = self.to_cycles()
        str = ""
        for e in tmp:
            # Если не перестановка вида (i)
            if e.__str__().__len__() != 3:
                str = str + e.__str__()
        return str


class Cycle:
    def __init__(self, data=None):
        # Для того, чтобы понимать границы
        self.max_el = 0
        if data is None:
            data = []
        self.data = [None]
        for el in data:
            if (el != "(") and (el != ")") and (el != " "):
                self.data.append(int(el))
                self.max_el = max(self.max_el, int(el))

    # Перевод цикла в перестановку
    def to_permutation(self):
        tmp = Permutation(self.max_el + 1)
        for i in range(1, len(self.data)):
            if i == len(self.data) - 1:
                tmp.perm[self.data[i]] = self.data[1]
            else:
                tmp.perm[self.data[i]] = self.data[i+1]
        return tmp

    def __str__(self):
        str = "("
        for e in self.data:
            if e != None:
                str = str + e.__str__() + ' '
        str_res = str[:-1]
        return str_res + ')'
