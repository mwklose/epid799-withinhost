from typing import Callable, List
import numpy as np
from src.person_objects.affinity_maturation_model import (
    AffinityMaturationModel,
)

from src.viral_objects.virus import Virus

# First and last year under study
FIRST_YEAR = 1968
LAST_YEAR = 2023


class Person:
    # have default constructor
    def __init__(
        self,
        id: str = "",
        birth_year: int = FIRST_YEAR,
        covariate_vector: List = [],
        infection_history: List = [],
    ):
        if len(covariate_vector) != len(infection_history):
            raise Exception("Covariates and Infections must be same length.")

        self.set_id(id)
        self.set_birth_year(birth_year)

        if not covariate_vector:
            self.set_covariate_by_function(Person.default_covariate_function)
        else:
            self.set_covariate_by_vector(covariate_vector=covariate_vector)

        if not infection_history:
            self.set_infection_history(Person.default_infection_history_function)
        else:
            self.set_infection_history(infection_history)
            self.infection_history_description = (
                f"Infection history by list: {infection_history}"
            )

        self.maturation_model = AffinityMaturationModel()

    def set_id(self, id: str) -> None:
        self.id = id

    def set_birth_year(self, birth_year: int) -> None:
        if birth_year < FIRST_YEAR or birth_year > LAST_YEAR:
            raise Exception(
                f"Invalid birth year; must be >= {FIRST_YEAR}, <= {LAST_YEAR}."
            )
        self.birth_year = birth_year

    def set_covariate_by_vector(self, covariate_vector: List) -> None:
        self.covariate = covariate_vector

    def set_covariate_by_function(
        self, covariate_function: Callable[[], int], *args
    ) -> None:
        # TODO: find solution where covariate function is single number versus more complex
        self.covariate = covariate_function(*args)

    def get_covariate_by_year(self, t: int) -> float:
        # TODO: add bounds checking
        return self.covariate[t]

    def set_infection_history(
        self,
        infection_history: List,
        force_of_infection_vector: List = [],
    ) -> None:
        # TODO: change to where infection history list is generated, rather than function calls
        if not force_of_infection_vector:
            self.infection_history = infection_history
        else:
            self.infection_history = infection_history_function(
                force_of_infection_vector
            )

    def get_infection_by_year(self, t: int) -> int:
        if t >= len(self.infection_history):
            return 0
        return self.infection_history[t]

    # Expose individual to the given year's virus.
    def expose_to_virus(self, virus: Virus) -> float:
        return self.maturation_model.exposure_to_virus(virus)

    def collect_time_to_clear(self) -> List:
        return self.maturation_model.collect_time_to_clear_infection()

    def collect_memory_cells_by_year(self) -> List:
        return self.maturation_model.collect_total_memory_cells_by_year()

    # ---------
    # Default functions for class, which we can make more complex or with more built-ins for the future.
    def default_covariate_function(*args) -> List[float]:
        return [0] * (LAST_YEAR - FIRST_YEAR + 1)

    def default_infection_history_function(
        *args,
        force_of_infection_vector: List[int],
    ) -> List[int]:
        rng = np.random.default_rng()
        random_pulls = rng.binomial(
            1, force_of_infection_vector, len(force_of_infection_vector)
        )

        return random_pulls

    def __str__(self) -> str:
        return (
            f"[{self.id}] -> b.{self.birth_year}, {self.infection_history_description}"
        )
