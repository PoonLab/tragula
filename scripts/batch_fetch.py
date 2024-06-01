from Bio import Entrez
from fetch import fetch, extract
import argparse
from csv import DictReader
import sys
from time import sleep
import os
import json

description = """
Retrieve all authors from PubMed using entries under 'query' header in a 
tab-separated input file.  Pause between queries to avoid spamming API.
"""

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument(
        "infile", type=argparse.FileType('r'),
        help="Tab-separated file to process")
    parser.add_argument(
        "dest", type=str, help="Directory to write JSON files")
    parser.add_argument(
        "--email", type=str,
        help="E-mail address, required for NCBI API transactions")
    parser.add_argument(
        "--retmax", type=int, default=200,
        help="Limit number of articles to retrieve per author.")
    parser.add_argument(
        "--verbose", action="store_true", help="Print messages.")
    parser.add_argument(
        "--overwrite", action="store_true", help="Replace existing JSON files.")
    args = parser.parse_args()

    if not os.path.exists(args.dest):
        sys.stderr.write(f"Error: destination path {args.dest} does not exist!\n")
        sys.exit()

    reader = DictReader(args.infile, delimiter="\t")
    if "query" not in reader.fieldnames:
        sys.stderr.write("ERROR: input file does not contain 'query' in header\n")
        sys.exit()

    Entrez.email = args.email
    for row in reader:
        query = row["query"]
        if query == "":
            if args.verbose:
                sys.stderr.write("Skipping empty query...\n")
            sleep(1)
            continue

        # construct filename from author names
        lastname = row['lastname'].replace(' ', '_')
        firstname = row['forename'].replace(' ', '_')
        fn = f"{lastname}_{firstname}.json"
        outfile = os.path.join(args.dest, fn)
        if os.path.exists(outfile) and not args.overwrite:
            # avoid unnecessary calls to API
            sys.stderr.write(f"File {outfile} exists, use --overwrite to replace.\n")
            continue

        # retrieve articles
        articles = fetch(row["query"], retmax=args.retmax)
        records = [extract(art) for art in articles]
        records = [r for r in records if r is not None]
        if args.verbose:
            sys.stderr.write(f"Query \"{query}\" retrieved {len(records)} "
                             f"records (limit {args.retmax})\n")

        # write results to JSON file
        with open(outfile, 'w') as handle:
            json.dump(records, handle, indent=2)

        sleep(1)  # pause
