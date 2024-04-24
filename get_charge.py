#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys


def main(file_path):
    f = open(file_path, 'r')
    lines = f.readlines()
    for line in lines:
        words = line.split()
        print(words)
    # Return an integer exit status
    return result


if __name__ == "__main__":
    # if len(sys.argv) != 2:
    #     print("Usage: python my_script.py file_path")
    #     sys.exit(1)  # Return a non-zero exit status to indicate an error

    # # Get the file path from the command-line argument
    # file_path = sys.argv[1]

    # Call the function to process the file and return the result
    result = main('residues.txt')
