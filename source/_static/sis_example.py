import random
import pandas as pd
from plotnine import *

import jsf

random.seed(7)

x0 = [1000 - 2, 2]
rates = lambda x, _: [2e-3 * x[0] * x[1], 1.0 * x[1]]
reactant_matrix = [[1, 1], [0, 1]]
product_matrix = [[0, 2], [1, 0]]

t_max = 10.0


stoich = {
    "nu": [
        [a - b for a, b in zip(r1, r2)]
        for r1, r2 in zip(product_matrix, reactant_matrix)
    ],
    "DoDisc": [1, 1],
    "nuReactant": reactant_matrix,
    "nuProduct": product_matrix,
}
my_opts = {"EnforceDo": [0, 0], "dt": 0.1, "SwitchingThreshold": [50, 50]}


sim = jsf.jsf(x0, rates, stoich, t_max, config=my_opts, method="operator-splitting")

sim_df = pd.DataFrame(
    {"time": sim[1], "susceptible": sim[0][0], "infectious": sim[0][1]}
).melt(id_vars=["time"], value_vars=["susceptible", "infectious"])

sim_p9 = (
    ggplot()
    + geom_line(data=sim_df, mapping=aes(x="time", y="value", colour="variable"))
    + geom_hline(yintercept=my_opts["SwitchingThreshold"][1], linetype="dashed")
    + scale_y_sqrt(name="Population size")
    + labs(x="Time", colour="Status")
    + theme(legend_position="top")
    + theme_bw()
)

sim_p9.save("sis_example.png", height=4, width=6)
