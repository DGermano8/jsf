import pandas as pd
import random
import time
import matplotlib.pyplot as plt
import jsf.__initOriginal__ as __initOriginal__
import jsf.__init__ as __init__

import numpy as np

random.seed(42)

x0 = [50, 10]

mA = 2.00
mB = 0.05
mC = 1.50

rates = lambda x, _: [mA * x[0],
                      mC * x[1],
                      mB * x[0] * x[1]]

reactant_matrix = [[1 , 0],
                   [0 , 1],
                   [1 , 1]]

product_matrix = [[2 , 0],
                  [0 , 0],
                  [0 , 2]]

t_max = 10

stoich = {

"nu": [
        [a - b for a, b in zip(r1, r2)] for r1, r2 in zip(product_matrix, reactant_matrix)
      ],

    "DoDisc": [1, 1],
    "nuReactant": reactant_matrix,
    "nuProduct": product_matrix,
}

my_opts = {"EnforceDo": [0, 0], "dt": 0.001, "SwitchingThreshold": [30, 30]}

random.seed(42)
start = time.perf_counter()
sim_1 = __initOriginal__.jsf(x0, rates, stoich, t_max, config=my_opts, method="operator-splitting")
end = time.perf_counter()
print(f"sim_1 execution time: {end - start:.4f} seconds")

random.seed(42)
start = time.perf_counter()
sim_2 = __init__.jsf(x0, rates, stoich, t_max, config=my_opts, method="operator-splitting")
end = time.perf_counter()
print(f"sim_2 execution time: {end - start:.4f} seconds")


# plt.plot(sim_1[1], sim_1[0][0], label="Prey", color=(0.9294, 0.4824, 0.4824), linewidth=2.5, marker='.', markersize=0.0)
# plt.plot(sim_1[1], sim_1[0][1], label="Predator", color=(0.5137, 0.3765, 0.5882),     linewidth=2.5, marker='.', markersize=0.0)

# plt.plot(sim_2[1], sim_2[0][0], label="Prey", color=(0.9294, 0.4824, 0.4824), linewidth=2.5, linestyle='--', markersize=0.0)
# plt.plot(sim_2[1], sim_2[0][1], label="Predator", color=(0.5137, 0.3765, 0.5882),     linewidth=2.5, linestyle='--', markersize=0.0)

plt.plot(np.array(sim_1[1]), np.array(sim_1[0][0])-np.array(sim_2[0][0]), label="Prey", color=(0.9294, 0.4824, 0.4824), linewidth=2.5, marker='.', markersize=0.0)
plt.plot(np.array(sim_1[1]), np.array(sim_1[0][1])-np.array(sim_2[0][1]), label="Predator", color=(0.5137, 0.3765, 0.5882),     linewidth=2.5, marker='.', markersize=0.0)

# plt.plot(sim_2[1], sim_2[0][0], label="Prey", color=(0.9294, 0.4824, 0.4824), linewidth=2.5, linestyle='--', markersize=0.0)
# plt.plot(sim_2[1], sim_2[0][1], label="Predator", color=(0.5137, 0.3765, 0.5882),     linewidth=2.5, linestyle='--', markersize=0.0)



plt.axhline(y=my_opts["SwitchingThreshold"][1], color=(0.5, 0.5, 0.5), linestyle="--")
plt.xlabel("Time")
plt.ylabel("Population size")
plt.legend()
plt.savefig('PredatorPrey.png', dpi=500)
plt.show()

# save matplotlib plot to file
