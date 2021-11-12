import numpy as np
from scipy import stats

# Hypothesis
lifetime_hours = 1000  # at least
alpha = 0.02

# Sample data
n = 16
x_bar_hours = 987.5
var_hours = 400
std_hours = np.sqrt(var_hours)
print(f"{std_hours=}")

t_stat = (x_bar_hours-lifetime_hours) / (std_hours/np.sqrt(n))
df = n-1

# t_stat = 2.2361
# df = 19
p = stats.t.cdf(t_stat, df)

print(f"{t_stat=}")
print(f"{p=} area to left of t_stat")
print(f"{1-p=} area to right of t_stat")

cv = stats.t.ppf(1.0 - alpha, df)
print(f"{cv=}")

