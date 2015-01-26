import unittest

from unassign.util import (
    uniq,
    )


class FunctionTests(unittest.TestCase):
    def test_uniq(self):
        """uniq should preserve the order of the input list."""
        a = [5, 5, 4, 5, 3, 5, 2, 1, 2, 3]
        self.assertEquals(uniq(a), [5, 4, 3, 2, 1])


if __name__ == "__main__":
    unittest.main()

