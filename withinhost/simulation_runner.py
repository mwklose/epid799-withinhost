from typing import List, Callable, Tuple

from src.person_objects.population import Population
from src.viral_objects.virus import Virus

import matplotlib.pyplot as plt


class SimulationRunner:
    def __init__(self, file_string: str, drift_function: Callable = lambda t: 5 * t):
        self.population = Population(file_string=file_string)
        self.year_range, self.virus_properties = self.generate_virus_history(
            drift_function=drift_function
        )

    def generate_virus_history(self, drift_function: Callable) -> Tuple[range, List]:
        years_needed = max(
            [len(person.infection_history) for person in self.population.list_of_people]
        )

        viral_history = [
            Virus(100, drift_function(t=time)) for time in range(0, years_needed)
        ]

        return range(0, years_needed), viral_history

    def run_simulation(self):
        for year in self.year_range:
            self.population.expose_to_virus(year, self.virus_properties[year])

    def visualize_simulation(self, img_name):
        plt.style.use("bmh")
        time_to_clear = self.population.collect_time_to_clear()
        NMODELS = 5
        # Plot 1: time to clear event.
        fig, ax = plt.subplots()
        for strategy, ttc in enumerate(time_to_clear):
            year = [(strategy + 1) * i for i in range(0, NMODELS)] + [30]
            time = ttc.values()

            ax.plot(
                year, time, label=f"Strategy {strategy + 1}", linestyle="-", alpha=0.5
            )
            ax.scatter(year, time)
            ax.set_xlabel("Year")
            ax.set_ylabel("Time to clear infection")

        ax.legend(loc=4, ncols=2)
        ax.set_title("Time to clear infection, by Strategy")
        plt.show()

        # Plot 2: Memory Cell Count. Best displayed by stacked plot.
        NCOL = 2
        NROW = 3
        fig, ax2 = plt.subplot_mosaic(
            [
                ["Strategy 1", "Strategy 2"],
                ["Strategy 3", "Strategy 4"],
                ["Strategy 5", "."],
            ],
            sharey=True,
        )

        total_memory_cell_count = self.population.collect_memory_cells_by_year()
        # Iterate through dictionary for each strategy
        for strategy, memory_cell in enumerate(total_memory_cell_count):
            year = [(strategy + 1) * i for i in range(0, NMODELS)] + [30]
            # Create stacked plot data structure.
            # Initialize empty dictionaries within each strategy
            plot_dict = {
                f"Infection {model + 1}": [] for model, _ in enumerate(memory_cell)
            }
            # Go through each cell model and append it.
            for num, value in enumerate(memory_cell.values()):
                padded_values = value + ([0] * (NMODELS + 1 - len(value)))
                for i, m in enumerate(padded_values):
                    plot_dict[f"Infection {i + 1}"].append(m)

            this_plot = ax2[f"Strategy {strategy + 1}"].stackplot(
                year, plot_dict.values(), labels=plot_dict.keys()
            )

            ax2[f"Strategy {strategy + 1}"].set_title(f"Strategy {strategy + 1}")
        ax2["Strategy 5"].legend(
            handles=this_plot,
            loc="center right",
            ncol=3,
            bbox_to_anchor=(2, 0.5),
        )
        plt.show()


if __name__ == "__main__":
    sim = SimulationRunner("birth_data.csv")
    sim.run_simulation()
    sim.visualize_simulation("output.png")

    # Simulation 2; drift function as alternation between
    sim2 = SimulationRunner(
        "birth_data.csv", drift_function=lambda t: 100 * ((t % 8) >= 4) + 0.1 * t
    )
    sim2.run_simulation()
    sim2.visualize_simulation("output2.png")
