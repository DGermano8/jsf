import numpy as np
import matplotlib.pyplot as plt

def birth_death_simulation(birth_rate, death_rate, initial_population, time_steps):
    population = [initial_population]

    for t in range(1, time_steps):
        current = population[-1]

        births = np.random.poisson(birth_rate * current)
        deaths = np.random.poisson(death_rate * current)

        new_population = current + births - deaths

        if new_population < 0:
            new_population = 0

        population.append(new_population)

    return population
if __name__ == "__main__":
    growth = birth_death_simulation(0.1, 0.05, 50, 100)
    decline = birth_death_simulation(0.05, 0.1, 50, 100)

    print("Growth:", growth)
    print("Decline:", decline)
if __name__ == "__main__":
    growth = birth_death_simulation(0.1, 0.05, 50, 100)
    decline = birth_death_simulation(0.05, 0.1, 50, 100)

    plt.plot(growth, label="Birth > Death (Growth)")
    plt.plot(decline, label="Death > Birth (Decline)")

    plt.xlabel("Time")
    plt.ylabel("Population")
    plt.title("Birth-Death Process Simulation")
    plt.legend()

    plt.savefig("birth_death_plot.png")  # saves image
    plt.show()