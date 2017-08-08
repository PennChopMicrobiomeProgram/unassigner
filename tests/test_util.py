import unittest

from unassign.util import *

class FunctionTests(unittest.TestCase):
    def test_uniq(self):
        """uniq should preserve the order of the input list."""
        self.assertEqual(
            uniq([5, 5, 4, 5, 3, 5, 2, 1, 2, 3, 3]), [5, 4, 3, 2, 1])

    def test_count_while_equal(self):
        self.assertEqual(count_while_equal("----a--sd-fb", "-"), 4)

    def test_count_while_equal_nomatch(self):
        self.assertEqual(count_while_equal("----a--sd-fb", "a"), 0)

    def test_count_while_equal_reversed(self):
        self.assertEqual(count_while_equal(reversed("----asd-fbbb"), "b"), 3)

    def test_count_matching_pairs(self):
        self.assertEqual(
            count_matching_pairs(zip("ABCDEFG", "ABCDEFG")), (7, 7))
        self.assertEqual(
            count_matching_pairs(zip("ABCDEFG", "ABXDYFZ")), (4, 7))

    def test_process_exits_successfully(self):
        self.assertEqual(process_exits_successfully(["echo", "hello"]), True)

if __name__ == "__main__":
    unittest.main()
