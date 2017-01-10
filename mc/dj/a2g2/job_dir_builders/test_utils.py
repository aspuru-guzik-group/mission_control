import difflib
import filecmp
import os


def assert_dirs_equal(test=None, left=None, right=None):
    cmp_result = filecmp.dircmp(left, right)
    test.assertEqual(set(cmp_result.left_only), set())
    test.assertEqual(set(cmp_result.right_only), set())
    for differing_file in cmp_result.diff_files:
        fail_msg = "Files differ:\n%s" % diff(
            left_file=os.path.join(left, differing_file),
            right_file=os.path.join(right, differing_file))
        test.fail(fail_msg)

def diff(left_file=None, right_file=None):
    diff_lines = difflib.unified_diff(
        open(left_file).readlines(), open(right_file).readlines(),
        fromfile=left_file, tofile=right_file)
    return "".join([line for line in diff_lines])
