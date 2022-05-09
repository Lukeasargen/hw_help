import random

count = 3

# [chapter, problems in chapter]
problems = [
    [1, 63],
    [2, 53],
    [3, 141],
    [4, 74],
    [5, 114],
    [6, 67],
    [7, 116],
    [8, 102],
    [9, 92],
    [10, 58],
    [11, 73],
    [12, 108],
    [13, 102],
    # [14, 56],
]

for i in range(count):
    chp, m = random.choice(problems)
    p = random.randint(0,m)
    print(f"{chp}.{p}")
