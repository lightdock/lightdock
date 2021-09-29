#!/usr/bin/env python3

"""Translates a LightDock output file to CSV format"""

import argparse
import os
from lightdock.util.logger import LoggingManager


log = LoggingManager.get_logger("lightdock2csv")
DEFAULT_SEP = ","


def parse_command_line():
    parser = argparse.ArgumentParser(prog="lightdock2csv")
    parser.add_argument("ranking_file", help="lightdock ranking file name")
    parser.add_argument("csv_file", help="csv format file name")
    parser.add_argument(
        "-s",
        "--sep",
        help="csv file separator character",
        dest="csv_sep",
        default=DEFAULT_SEP,
    )
    args = parser.parse_args()
    return args


if __name__ == "__main__":
    try:
        # Parse command line
        args = parse_command_line()
        fin = open(args.ranking_file)
        lines = fin.readlines()
        fin.close()
        output = open(args.csv_file, "w")
        output.write(args.csv_sep.join(lines[0].split()) + os.linesep)
        count = 0
        for line in lines[1:]:
            first = line.index("(")
            last = line.index(")")
            fields = []
            fields.extend(line[:first].split())
            fields.append('"' + line[first : last + 1].strip() + '"')
            fields.extend(line[last + 1 :].split())
            output.write(args.csv_sep.join(fields) + os.linesep)
            count += 1
        output.close()
        log.info("%d lines written" % count)
        log.info("Done.")

    except KeyboardInterrupt:
        log.info("Caught interrupt...")
        log.info("bye.")
