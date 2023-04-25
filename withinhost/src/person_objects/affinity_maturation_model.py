from typing import List, Tuple, Callable, Dict
from src.viral_objects.virus import Virus

import numpy as np
from scipy.linalg import block_diag
from scipy.integrate import solve_ivp
import scipy
import itertools


# Virus doubles whenever possible
PLASMA_TO_MEMORY_FACTOR = 0.1
MEMORY_TO_PLASMA_FACTOR = 1
MEMORY_DECAY = 0
PLASMA_DECAY = 0.5


class AffinityMaturationModel:
    def __init__(self):
        self.antibody_models: List[AntibodyModel] = []
        self.time_to_clear: Dict = {}
        self.total_memory_cell_count: Dict = {}

    def exposure_to_virus(self, virus: Virus) -> float:
        """
        This method simulates a person's exposure to a virus and computes the person's immunity response.

        Args:
            virus (Virus): A Virus object that represents the virus the person is exposed to.

        Returns:
            exposure_results (float): A float that represents the person's immunity response to the virus.
        """
        # Case 1: person is actually not exposed.

        # Case 2: person is exposed.
        # 1. Add another model to the system for the current virus
        self.antibody_models.append(AntibodyModel(virus))

        # Ultimately, working towards setting up ODE to solve.
        # A. Collect baseline conditions
        # B. Collect differential equations for each model
        starting_values, differential_equations = self.construct_differential_equations(
            virus
        )
        init_values = list(starting_values)
        print(f"starting_values = {init_values}")

        # D. Set up ODE solver
        # E. Solve
        def virus_zero_cross(t, y):
            return y[0]

        def plasma_zero_cross(t, y):
            return y[-2]

        virus_zero_cross.direction = -1
        plasma_zero_cross.terminal = True
        plasma_zero_cross.direction = -1

        ode_solution = solve_ivp(
            fun=differential_equations,
            t_span=[0, 100],
            y0=init_values,
            events=[virus_zero_cross, plasma_zero_cross],
        )
        # F. Get results
        print(ode_solution)
        # G. Write memory cell back to each model
        exposure_results = self.extract_ode_solution(ode_solution, virus)

        return exposure_results

    def extract_ode_solution(self, ode_soln, virus: Virus) -> float:
        """
        This method extracts the solution from the ODE solver and updates the memory cells for the given virus.

        Args:
            ode_soln: A solution object returned by the ODE solver.
            virus (Virus): A Virus object that represents the virus for which memory cells need to be updated.

        Returns:
            memory_cell_count (float): A float that represents the total number of memory cells for the given virus.
        """

        self.time_to_clear[virus.genetic_code] = ode_soln.t_events[0][0]
        self.total_memory_cell_count[virus.genetic_code] = self.write_memory_cells(
            ode_soln.y_events[1][0]
        )
        return sum(self.total_memory_cell_count[virus.genetic_code])

    def write_memory_cells(self, solution_array: List) -> None:
        """
        This method updates the memory cell counts for each antibody model using the solution array from the ODE solver.

        Args:
            solution_array (List): A list of floats that represents the solution array from the ODE solver.

        Returns:
            None
        """
        b_and_m_values = solution_array[1:]
        b_and_m_columns = [
            b_and_m_values[i : i + 2] for i in range(0, len(b_and_m_values), 2)
        ]
        model_columns = zip(self.antibody_models, b_and_m_columns)

        for model, (b, m) in model_columns:
            model.b_cells = 0
            model.m_cells = m

        return [model.m_cells for model in self.antibody_models]

    def construct_differential_equations(self, virus: Virus) -> Callable:
        # Initial parameters, flatten to 1D list/array

        initial_parameters = [virus.viral_load]
        for model in self.antibody_models:
            initial_parameters = itertools.chain(
                initial_parameters, model.get_conditions(virus)
            )

        # Constructing differential equations themselves
        # t is required by solve_ivp, y is initialized with y0.
        def differential_equations(t, y):
            """
            Computes the differential equations for a viral response model.

            Parameters:
            t (float): The current time step.
            y (List[float]): The current values of the state variables.

            Returns:
            List[float]: The rates of change for each of the state variables.

            Notes:
            This function computes the rates of change for a viral response model
            consisting of a viral load and multiple types of B-cells (memory and plasma).
            The model parameters are stored in the state variable `y`.
            The differential equations are solved using a numerical integration method.

            """
            # Viral response - double current virus level,
            viral_duplication = max(0, y[0])

            model_params = [x for x in y[1:]]
            model_params = [
                model_params[i : i + 2] for i in range(0, len(model_params), 2)
            ]

            # Flatten to 1D array
            viral_response = [
                model.viral_response(current_viral_load=y[0], virus=virus, b_cells=b)
                for model, (b, m) in zip(self.antibody_models, model_params)
            ]

            viral_de = viral_duplication + sum(itertools.chain(viral_response))

            # B-cell responses: dB = V/D + XMV/D - decay
            # First, need to get model parameters, then can create pairs
            # Model needs viral load, memory cells, genetic code
            b_response = [
                model.plasma_response(
                    current_viral_load=y[0], virus=virus, b_cells=b, m_cells=m
                )
                for model, (b, m) in zip(self.antibody_models, model_params)
            ]

            m_response = [
                model.memory_response(
                    current_viral_load=y[0], virus=virus, b_cells=b, m_cells=m
                )
                for model, (b, m) in zip(self.antibody_models, model_params)
            ]

            # Make differential equations into 1D array
            zipped_responses = zip(b_response, m_response)
            responses = list(itertools.chain.from_iterable(zipped_responses))
            values = list(itertools.chain([viral_de], responses))
            return values

        return initial_parameters, differential_equations

    def save_model_params() -> List[Tuple[int, int]]:
        """
        Returns a list of tuples containing the virus genetic code and number of memory cells for each
            antibody model in the simulation.

        Returns:
            List of tuples containing the virus genetic code and number of memory cells for each
            antibody model in the simulation.
        """
        model_params = [
            (model.virus_genetic_code, model.m_cells) for model in self.antibody_models
        ]
        return model_params

    def collect_time_to_clear_infection(self):
        return self.time_to_clear

    def collect_total_memory_cells_by_year(self):
        return self.total_memory_cell_count


# Create another class just to hold the components for a certain virus in a certain year.
class AntibodyModel:
    def __init__(self, virus: Virus):
        self.virus_genetic_code = virus.genetic_code
        self.b_cells = 0
        self.m_cells = 0

    def get_conditions(self, virus: Virus) -> List[float]:
        """
        Returns the conditions of a Virus for a given antibody model.

        Args:
        - virus: A Virus object representing the virus.

        Returns:
        - A list of floats containing the number of B-cells and memory cells for the given Virus.
        """

        if self.virus_genetic_code == virus.genetic_code:
            b = 0
        else:
            b = self.b_cells

        m = self.m_cells

        return [b, m]

    def viral_response(
        self, current_viral_load: float, virus: Virus, b_cells: float
    ) -> List[float]:
        """
        This function calculates the viral response to a given virus, based on the current viral load and number of B-cells. It returns a list of floats representing the viral response.

        Args:
        - current_viral_load (float): The current viral load
        - virus (Virus): The virus to which the response is being calculated
        - b_cells (float): The number of B-cells

        Returns:
        - List[float]: A list of floats representing the viral response. The list contains only one element, which is the viral response based on the current viral load and number of B-cells. The value is calculated as follows:
        """

        b_response = (
            -1 * b_cells / (1 + virus.get_genetic_distance(self.virus_genetic_code))
        )
        # Keep flexibility for later; memory cell response can be non-zero.
        m_response = 0
        return (current_viral_load > 0) * b_response + m_response

    def plasma_response(
        self, current_viral_load: float, virus: Virus, b_cells: float, m_cells: float
    ) -> float:
        """
        Function Name: plasma_response

        Input:

            current_viral_load (float): current level of viral load
            virus (Virus): instance of a Virus object representing the current virus
            b_cells (float): number of B cells
            m_cells (float): number of memory cells

        Output:

            dB (float): the plasma response value

        Functionality:
        This function calculates the plasma response of an antibody model, given the current viral load, virus, B cells, and memory cells. It first calculates the genetic distance between the current virus and the virus the antibody model was trained on. It then computes the plasma response, dB, as follows:

            The first term is the viral load divided by (1 + genetic distance), or the reduction of viral load due to antibodies.
            The second term is the plasma-to-memory factor multiplied by the number of B cells, which represents the decay of plasma cells.
            The third term is the memory-to-plasma factor multiplied by the number of memory cells and the current viral load divided by (1 + genetic distance), or the production of plasma cells from memory cells.
            The last term is the decay rate of plasma cells.

        Finally, the function returns the computed dB value.
        """
        distance = virus.get_genetic_distance(self.virus_genetic_code)
        dB = (
            max(0, current_viral_load) / (1 + distance)
            - PLASMA_TO_MEMORY_FACTOR * max(0, b_cells)
            + MEMORY_TO_PLASMA_FACTOR
            * m_cells
            * max(0, current_viral_load)
            / (1 + distance)
            - PLASMA_DECAY
        )
        return dB

    def memory_response(
        self, current_viral_load: float, virus: Virus, b_cells: float, m_cells: float
    ) -> float:
        """
        Calculates the change in memory cells in response to a virus.

        Parameters:
        -----------
        current_viral_load : float
            The current viral load.
        virus : Virus
            The virus being responded to.
        b_cells : float
            The number of B-cells.
        m_cells : float
            The number of memory cells.

        Returns:
        --------
        float
            The change in memory cells in response to the virus.
        """
        distance = virus.get_genetic_distance(self.virus_genetic_code)
        dM = (
            PLASMA_TO_MEMORY_FACTOR * max(0, b_cells)
            - MEMORY_TO_PLASMA_FACTOR
            * m_cells
            * max(0, current_viral_load)
            / (1 + distance)
            - MEMORY_DECAY
        )
        return dM
