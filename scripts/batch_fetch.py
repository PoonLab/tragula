from Bio import Entrez
from fetch import fetch, extract
import argparse
from csv import DictReader
import sys
from time import sleep

description = """
Retrieve all authors in a tab-separated file using entries under 'query' 
header from PubMed database
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("infile", type=argparse.FileType('r'),
                        help="Tab-separated file to process")
    parser.add_argument("--email", type=str,
                        help="E-mail address, required for PubMed API transactions")
    args = parser.parse_args()

    reader = DictReader(args.infile, delimiter="\t")
    if "query" not in reader.fieldnames:
        sys.stderr.write("ERROR: input file does not contain 'query' in header\n")
        sys.exit()
    Entrez.email = args.email
    for row in reader:
        