import unittest
import numpy as np
import pandas as pd
from withinhost.src.person_objects.population import Population
from withinhost.src.person_objects.person import Person

from withinhost.src.viral_objects.virus import Virus


class PopulationUnitTest(unittest.TestCase):
    def test_init_from_string(self):
        pop = Population("birth_data.csv")

        for item in pop.list_of_people:
            print(item)
            self.assertIsInstance(item, Person)
