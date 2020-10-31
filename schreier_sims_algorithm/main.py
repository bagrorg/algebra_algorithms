from permutation import *
from schreier_tree import *
from stabilizer_chain import *
from tests import *

# ----------------------------------------------ЗАДАНИЕ------------------------------------------------------
# Создание массивов для перестановок при делении колоды на 3 части
# 2 1 3
task_swap_f_s = []
for i in range(1, 19):
    task_swap_f_s.append([i, i + 18])

# 3 2 1
task_swap_f_t = []
for i in range(1, 19):
    task_swap_f_t.append([i, i + 36])

# 3 1 2
task_swap_t_on_top = []
for i in range(1, 19):
    task_swap_t_on_top.append([i, i + 18, i + 36])

# 2 3 1 - квадрат прошлой
task_swap_s_on_top = []

# 1 3 2
task_swap_s_t = []
for i in range(19, 37):
    task_swap_s_t.append([i, i + 18])

task_shtirlez = [[1, 54], [2, 53]]

# -------------------------------
# Конвертация массивов в перестановки

perm213 = Permutation(55)
for e in task_swap_f_s:
    perm213 = perm213 * Cycle(e).to_permutation()

perm321 = Permutation(55)
for e in task_swap_f_t:
    perm321 = perm321 * Cycle(e).to_permutation()

perm312 = Permutation(55)
for e in task_swap_t_on_top:
    perm312 = perm312 * Cycle(e).to_permutation()

perm231 = perm312 ** 2
perm132 = Permutation(55)
for e in task_swap_s_t:
    perm132 = perm132 * Cycle(e).to_permutation()

# Перестановка с двумя картами
perm_sht_arr = []
t = 54
while t > 0:
    perm_sht_arr.append(t)
    t -= 2
    if t == 0:
        t = 53
perm_sht = Cycle(perm_sht_arr).to_permutation()

# -------------------------------
# Создание множества на основе перестановок
ans = set()
ans.add(perm132)
ans.add(perm312)
ans.add(perm321)
ans.add(perm213)
ans.add(perm231)
ans.add(perm_sht)
ans_set = StabilizerChain(ans)

# Посмотрим на порядок группы
print(ans_set.group_ord())

# Перестановка, которую мы ищем
task = Permutation(55)
task.perm = [x for x in range(54, 0, -1)]

# Вывод ответа
print(ans_set.contains(Cycle(task).to_permutation()))
