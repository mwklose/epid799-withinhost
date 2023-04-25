import unittest

from withinhost.src.person_objects.person import Person


class PersonUnitTest(unittest.TestCase):
    def test_default_infection_history_function(self):
        result = Person.default_infection_history_function(
            force_of_infection_vector=[1, 0, 1]
        )
        self.assertTrue(
            all(x == y for x, y in zip(result, [1, 0, 1])),
            "Certainty of infection history",
        )

        long_result = Person.default_infection_history_function(
            force_of_infection_vector=[0.5] * 100
        )
        self.assertTrue(
            any(x == 1 for x in long_result), "Binomial choice works as expected"
        )

    def test_covariate_default(self):
        test_person = Person()

        self.assertEqual(
            0,
            test_person.get_covariate_by_year(0),
            "On construction, covariate is none",
        )

        test_person.set_covariate_by_function(Person.default_covariate_function)

        self.assertEqual(
            0,
            test_person.get_covariate_by_year(0),
            "Test covariate equal after setting with default",
        )
