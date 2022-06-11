import os


def compare_two_files(setup_file, expected_setup_file, ignore=None):
    """Compares two files and maybe using an ignore list of words"""
    lines = []
    with open(setup_file) as ih:
        for line in ih.readlines():
            if type(ignore) == list and any(x in line for x in ignore):
                continue
            else:
                lines.append(line.rstrip(os.linesep))

    expected_lines = []
    with open(expected_setup_file) as ih:
        for line in ih.readlines():
            if type(ignore) == list and any(x in line for x in ignore):
                continue
            else:
                expected_lines.append(line.rstrip(os.linesep))

    if len(lines) != len(expected_lines):
        return False

    for line1, line2 in zip(lines, expected_lines):
        if line1 != line2:
            return False

    return True
