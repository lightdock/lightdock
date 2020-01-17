#!/usr/bin/env python3

"""Calculates the success rate (Top 10 and Top100 of NN) for PyDock or LightDock results"""

import argparse
from lightdock.util.table import Table


DEFAULT_RMSD_CUTOFF = 10.0


class PydockData(object):
    def __init__(self, top10, top100, hits, min_rmsd_10, min_rmsd_100, min_rmsd, first_nn):
        self.top10 = top10
        self.top100 = top100
        self.hits = hits
        self.min_rmsd_10 = min_rmsd_10
        self.min_rmsd_100 = min_rmsd_100
        self.min_rmsd = min_rmsd
        self.first_nn = first_nn

    def show(self, csv=False):
        if csv:
            print("%d:%d:%d:%5.3f:%5.3f:%5.3f:%d" % (self.top10, self.top100, self.hits, self.min_rmsd_10,
                                                     self.min_rmsd_100, self.min_rmsd, self.first_nn))
        else:
            print("PyDock Top10: %d" % self.top10)
            print("PyDock Top100: %d" % self.top100)
            print("PyDock Total Hits: %d" % self.hits)
            print("PyDock Min(RMSD) Top 10: %5.3f" % self.min_rmsd_10)
            print("PyDock Min(RMSD) Top 100: %5.3f" % self.min_rmsd_100)
            print("PyDock Min(RMSD): %5.3f" % self.min_rmsd)
            print("PyDock First Near Native: %d" % self.first_nn)


class LightdockData(object):
    def __init__(self, top10, top100, hits, min_rmsd_10, min_rmsd_100, min_rmsd,
                 mean_clashes, min_clashes, max_clashes, mean_scoring_10, mean_scoring_100, solutions, first_nn):
        self.top10 = top10
        self.top100 = top100
        self.hits = hits
        self.min_rmsd_10 = min_rmsd_10
        self.min_rmsd_100 = min_rmsd_100
        self.min_rmsd = min_rmsd
        self.mean_clashes = mean_clashes
        self.min_clashes = min_clashes
        self.max_clashes = max_clashes
        self.mean_scoring_10 = mean_scoring_10
        self.mean_scoring_100 = mean_scoring_100
        self.solutions = solutions
        self.first_nn = first_nn

    def show(self, csv=False):
        if csv:
            print("%d:%d:%d:%5.3f:%5.3f:%5.3f:%5.3f:%5.3f:%5.3f:%5.3f:%5.3f:%d:%d" % (self.top10, self.top100,
                                                                                      self.hits, self.min_rmsd_10,
                                                                                      self.min_rmsd_100, self.min_rmsd,
                                                                                      self.mean_clashes,
                                                                                      self.min_clashes,
                                                                                      self.max_clashes,
                                                                                      self.mean_scoring_10,
                                                                                      self.mean_scoring_100,
                                                                                      self.solutions, self.first_nn))
        else:
            print("Lightdock Top 10: %d" % self.top10)
            print("Lightdock Top 100: %d" % self.top100)
            print("Lightdock Total Hits: %d" % self.hits)
            print("Lightdock Min(RMSD) Top 10: %5.3f" % self.min_rmsd_10)
            print("Lightdock Min(RMSD) Top 100: %5.3f" % self.min_rmsd_100)
            print("Lightdock Min(RMSD): %5.3f" % self.min_rmsd)
            print("Lightdock Mean Clashes: %5.3f" % self.mean_clashes)
            print("Lightdock Min Clashes: %5.3f" % self.min_clashes)
            print("Lightdock Max Clashes: %5.3f" % self.max_clashes)
            print("Lightdock Mean Scoring Top 10: %5.3f" % self.mean_scoring_10)
            print("Lightdock Mean Scoring Top 100: %5.3f" % self.mean_scoring_100)
            print("Lightdock number of solutions: %d" % self.solutions)
            print("First Near Native position: %d" % self.first_nn)


def parse_arguments():
    arg_parser = argparse.ArgumentParser(prog='hits')
    arg_parser.add_argument("energy_file", help="energy file to analyze",
                            type=str, metavar="energy_file")
    arg_parser.add_argument("-r", "--rmsd", help="RMSD to be considered as a hit",
                            dest="rmsd_cutoff", type=float, default=DEFAULT_RMSD_CUTOFF)
    arg_parser.add_argument("-p", "--pydock", help="energy file is considered as a pyDock output",
                            dest="pydock", action='store_true')
    arg_parser.add_argument("--csv", help="output result in csv format",
                            dest="csv", action='store_true')
    args = arg_parser.parse_args()
    return args


def read_pydock_data(energy_file, rmsd_cutoff):
    pydock = Table.read(energy_file)
    pydock_top10 = 0
    pydock_top100 = 0
    pydock_min_rmsd_10 = 999.99
    pydock_min_rmsd_100 = 999.99
    pydock_min_rmsd = 999.99
    pydock_hits = 0
    first_nn = 99999

    for i, rmsd in enumerate(pydock['RMSD']):
        rmsd = float(rmsd)
        if i < 10:
            if rmsd <= rmsd_cutoff:
                pydock_top10 += 1
            if rmsd < pydock_min_rmsd_10:
                pydock_min_rmsd_10 = rmsd
        if i < 100:
            if rmsd <= rmsd_cutoff:
                pydock_top100 += 1
            if rmsd < pydock_min_rmsd_100:
                pydock_min_rmsd_100 = rmsd
        if rmsd <= rmsd_cutoff:
            pydock_hits += 1
        if rmsd < pydock_min_rmsd:
            pydock_min_rmsd = rmsd
        if first_nn == 99999 and rmsd <= rmsd_cutoff:
            first_nn = i + 1

    pydock_data = PydockData(pydock_top10, pydock_top100, pydock_hits, pydock_min_rmsd_10,
                             pydock_min_rmsd_100, pydock_min_rmsd, first_nn)
    return pydock_data


def read_lightdock_data(energy_file, rmsd_cutoff):
    lightdock_file = open(energy_file)
    lines = lightdock_file.readlines()[1:]
    lightdock_file.close()

    top10 = 0
    top100 = 0
    hits = 0
    mean_clashes = 0.0
    max_clashes = -20000
    min_clashes = 20000
    min_rmsd_10 = 999.9
    min_rmsd_100 = 999.9
    min_rmsd = 999.9
    mean_scoring_10 = 0.0
    mean_scoring_100 = 0.0
    first_nn = len(lines)
    for i, line in enumerate(lines):
        fields = line.rstrip().split()
        rmsd = float(fields[-4])
        clashes = float(fields[-2])
        scoring = float(fields[-1])
        if rmsd <= rmsd_cutoff and first_nn == len(lines):
            first_nn = i + 1 
        if i < 10:
            if rmsd < min_rmsd_10:
                min_rmsd_10 = rmsd
            if rmsd <= rmsd_cutoff:
                top10 += 1
            mean_scoring_10 += scoring

        if i < 100:
            if rmsd <= rmsd_cutoff:
                top100 += 1
            if rmsd < min_rmsd_100:
                min_rmsd_100 = rmsd
            mean_scoring_100 += scoring

        if rmsd <= rmsd_cutoff:
            hits += 1

        mean_clashes += clashes
        if clashes < min_clashes:
            min_clashes = clashes
        if clashes > max_clashes:
            max_clashes = clashes

        if rmsd < min_rmsd:
            min_rmsd = rmsd

    mean_clashes /= len(lines)
    mean_scoring_10 /= 10.0
    mean_scoring_100 /= 100.0
    lightdock_data = LightdockData(top10, top100, hits, min_rmsd_10, min_rmsd_100, min_rmsd,
                                   mean_clashes, min_clashes, max_clashes,
                                   mean_scoring_10, mean_scoring_100, len(lines), first_nn)
    return lightdock_data


if __name__ == "__main__":

    parser = parse_arguments()

    if parser.pydock:
        data = read_pydock_data(parser.energy_file, parser.rmsd_cutoff)
    else:
        data = read_lightdock_data(parser.energy_file, parser.rmsd_cutoff)

    data.show(parser.csv)
