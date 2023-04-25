import unittest
import numpy as np
from withinhost.src.person_objects.affinity_maturation_model import (
    AffinityMaturationModel,
    AntibodyModel,
)

from withinhost.src.viral_objects.virus import Virus


class AntibodyModelUnitTest(unittest.TestCase):
    def test_get_conditions(self):
        am = AntibodyModel(Virus(150, 5.0))

        other_virus = Virus(200, 15)

        self.assertTrue(all(x == 0 for x in am.get_conditions(other_virus)))

    def test_viral_response(self):
        virus1 = Virus(200, 15)
        virus2 = Virus(100, 20)

        am = AntibodyModel(virus1)
        am.b_cells = 100

        total = am.viral_response(virus2.viral_load, virus2, am.b_cells)

        self.assertEqual(total, -100 / 6)

    def test_plasma_response(self):
        virus1 = Virus(300, 10)
        virus2 = Virus(200, 30)

        am = AntibodyModel(virus1)
        plasma = am.plasma_response(virus2.viral_load, virus2, 100, 50)

        self.assertEqual(plasma, 200 * (1 + 1 * 50) / 21 - 0.5 - 0.1 * 100)

    def test_memory_response(self):
        virus1 = Virus(300, 10)
        virus2 = Virus(200, 30)

        am = AntibodyModel(virus1)
        am.b_cells = 100
        memory = am.memory_response(virus2.viral_load, virus2, 100, 50)

        self.assertEqual(memory, 0.1 * 100 - (50 * 200) / 21)

        memory = am.memory_response(virus1.viral_load, virus1, 100, 50)

        self.assertEqual(memory, 0.1 * 100 - (50 * 300))


class AffinityMaturationModelUnitTest(unittest.TestCase):
    def test_exposure(self):
        amm = AffinityMaturationModel()
        virus1 = Virus(100, 10)
        virus2 = Virus(100, 20)
        virus3 = Virus(100, 21)
        amm.exposure_to_virus(virus1)
        amm.exposure_to_virus(virus2)
        amm.exposure_to_virus(virus3)
