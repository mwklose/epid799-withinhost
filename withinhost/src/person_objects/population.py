import pandas as pd
from src.viral_objects.virus import Virus
from src.person_objects.person import Person


class Population:
    # TODO next:
    # since we are comparing treatment strategies, each "person" represents a treatment strategy
    def __init__(self, file_string: str):
        """
        Expects data to be read in as a long format.
        """
        self.list_of_people = []
        df = pd.read_csv(file_string)

        for name, info in df.groupby(by="ID"):
            birth_year = info["Year"].min()
            covariate_vector = info["Covariate"].to_list()
            infection_vector = info["Infection"].to_list()
            self.list_of_people.append(
                Person(
                    id=name,
                    birth_year=birth_year,
                    covariate_vector=covariate_vector,
                    infection_history=infection_vector,
                )
            )

    def expose_to_virus(self, t: int, virus: Virus):
        # TODO: vectorization/speedup probably possible.
        for person in self.list_of_people:
            if person.get_infection_by_year(t):
                person.expose_to_virus(virus)

    def collect_time_to_clear(self):
        return [person.collect_time_to_clear() for person in self.list_of_people]

    def collect_memory_cells_by_year(self):
        return [person.collect_memory_cells_by_year() for person in self.list_of_people]

    def make_default_file(file_destination: str):
        pass

    def generate_random_population(*args):
        # TODO: generate random population with certain covariate distributions.
        pass
