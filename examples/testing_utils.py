import difflib
import filecmp
import glob
import json
import os
import re


def get_dir_tree_str(dir=None, spacer='-'):
    tree_str = ''
    for dirpath, dirnames, filenames in os.walk(dir):
        path = dirpath.split(os.sep)
        tree_str += "{spacing} {path}\n".format(
            spacing=((len(path) - 1) * spacer),
            path=os.path.basename(dirpath)
        )
        for filename in filenames:
            tree_str += "{spacing} {filename}\n".format(
                spacing=(len(path) * spacer),
                filename=filename
            )
    return tree_str


def get_dir_items(dir):
    items = glob.glob(dir + '/**', recursive=True)
    dir_prefix_len = len(dir) + 1
    relative_items = [item[dir_prefix_len:] for item in items]
    return sorted(relative_items[1:])


def assert_dirs_equal(test=None, left=None, right=None, ignore_patterns=None,
                      json_patterns=None):
    ignore_patterns = ignore_patterns or []
    json_patterns = json_patterns or [r'\.json$']
    cmp_result = filecmp.dircmp(left, right)
    test.assertEqual(set(cmp_result.left_only), set())
    test.assertEqual(set(cmp_result.right_only), set())
    for differing_file in cmp_result.diff_files:
        for ignore_pattern in ignore_patterns:
            if re.search(ignore_pattern, differing_file):
                continue
        left_file = os.path.join(left, differing_file)
        right_file = os.path.join(right, differing_file)
        check_as_json = False
        for json_pattern in json_patterns:
            if re.search(json_pattern, differing_file):
                check_as_json = True
                break
        if check_as_json:
            assert_json_files_equal(test=test, left_file=left_file,
                                    right_file=right_file)
        else:
            test.fail("Files differ:\n%s" % (
                diff(left_file=left_file, right_file=right_file)))


def assert_json_files_equal(test=None, left_file=None, right_file=None):
    test.assertEqual(json.load(open(left_file)), json.load(open(right_file)))


def diff(left_file=None, right_file=None):
    diff_lines = difflib.unified_diff(
        open(left_file).readlines(), open(right_file).readlines(),
        fromfile=left_file, tofile=right_file)
    return "".join([line for line in diff_lines])
