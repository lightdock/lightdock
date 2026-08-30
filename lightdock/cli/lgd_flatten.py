"""Transforms a saved ANM NumPy matrix from 3D to 1D"""

import numpy as np
import sys


def main():
    if len(sys.argv[1:]) != 2:
        print(f"Usage: {sys.argv[0]} input.npy output.npy")
        raise SystemExit("Wrong command line")

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    n = np.load(input_file)
    print("{} -> 1d".format(n.shape))
    np.save(output_file, n.flatten())


if __name__ == "__main__":
    main()
