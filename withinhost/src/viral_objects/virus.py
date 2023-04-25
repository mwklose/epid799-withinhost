from typing import List


class Virus:
    """
    The code defines two classes, `Virus` and `NullVirus`.

    ### `Virus` Class

    The `Virus` class has the following attributes and methods:

    #### Attributes:

    - `viral_load`: an integer representing the amount of virus in the body.
    - `genetic_code`: a float representing the genetic code of the virus.

    #### Methods:

    - `__init__(self, viral_load: int, genetic_code: float)`: The constructor method which initializes `viral_load` and `genetic_code` attributes.
    - `get_genetic_distance(self, genetic_code: float) -> float`: A method that calculates the genetic distance between the virus instance and another virus instance using their `genetic_code` attributes.
    - `__str__(self) -> str`: A method that returns a string representation of the `Virus` instance. It returns a string containing the `viral_load` and `genetic_code` attributes.
    """

    def __init__(self, viral_load: int, genetic_code: float):
        self.viral_load = viral_load
        self.genetic_code = genetic_code

    def get_genetic_distance(self, genetic_code: float) -> float:
        return abs(self.genetic_code - genetic_code)

    def __str__(self) -> str:
        return f"VL={self.viral_load} GC={self.genetic_code}"

    def generate_virus_from_file(file_path: str) -> List:
        pass
