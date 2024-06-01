import argparse
import json
import nltk
import os
from glob import glob
import sys
from nltk.stem.wordnet import WordNetLemmatizer
import re
import csv

description = """
Tokenize article text for all author-specific JSON files in a given 
folder.  Export word co-occurrence as sparse matrix, and author 
word frequency vectors as a second file.  If --counts and --matrix are 
not specified, then the script will write all word frequencies to the path
specified in --index (defaults to stdout).
"""

# tags to exclude
# https://www.ibm.com/docs/en/wca/3.5.0?topic=analytics-part-speech-tag-sets
xtags = {
    "IN",  # preposition / subordinating conjunction, e.g., "because"
    "CD",  # cardinal number, e.g., "1"
    "DT",  # determiner, e.g., "the", "my"
    "CC",  # coordinating conjunction, e.g., "and"
    "TO",  # preposition, "to""
    ":", "(", ")", ",", "."  # punctuation
}

wordnet = WordNetLemmatizer()  # convert plural to singular
numeric = re.compile("^-?[0-9]+\.?[0-9]*$")


def load_xwords(path):
    """ Import a set of words to exclude from file """
    xwords = set()
    with open(path) as handle:
        for line in handle:
            xwords.add(line.strip())
    return xwords


def process(files, xwords, debug=False, min_records=3):
    """
    Enumerate word co-occurrence across all texts, and track
    author-specific word frequencies.
    :param files:  list, paths to JSON files containing
    :param xwords:  set, words to exclude
    :param debug:  bool, get global word counts only
    :param min_count:  int, minimum number of records to proceed with analysis
    :returns:
      - by_author - dict, word counts by author
      - cooccur - dict, stores numbers of texts in which both words appear
    """
    by_author = {}
    by_document = []
    all_words = {}
    for f in files:
        author = os.path.basename(f).split('.')[0]
        if author in by_author:
            sys.stderr.write(f"ERROR: duplicate author {author} detected!")
            sys.exit()
        by_author.update({author: {}})

        with open(f) as handle:
            records = json.load(handle)

        if len(records) < min_count:
            sys.stderr.write(f"Fewer than {len(records}} records, skipping.\n")
            continue

        for r in records:
            text = f"{r['title']} {r['abstract']} {' '.join(r['keywords'])}"
            tokens = nltk.word_tokenize(text)
            tagged = nltk.pos_tag(tokens)  # part-of-speech tagger
            filtered = []
            for word, tag in tagged:
                if tag in xtags:
                    continue
                word = wordnet.lemmatize(word.lower())
                word = word.strip("'")
                if word in xwords or numeric.findall(word):
                    continue
                filtered.append(word)

            # increment word counts
            for word in filtered:
                # global counts
                if word not in all_words:
                    all_words.update({word: 0})
                all_words[word] += 1

                if debug:
                    continue  # skip next steps
                # count by author
                if word not in by_author[author]:
                    by_author[author].update({word: 0})
                by_author[author][word] += 1

            # record co-occurrences
            if not debug:
                by_document.append(
                    dict([(w, filtered.count(w)) for w in set(filtered)])
                )

    return all_words, by_author, by_document


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description)
    parser.add_argument(
        "indir", type=str, help="Path to folder containing JSON files.")
    parser.add_argument(
        "--counts", type=argparse.FileType('w'), required=False,
        help="Path to write by-author word counts (JSON)")
    parser.add_argument(
        "--matrix", type=argparse.FileType('w'), required=False,
        help="Path to write co-occurrence sparse matrix (CSV)")
    parser.add_argument(
        "--index", type=argparse.FileType('w'), default=sys.stdout,
        help="Path to write indexed words for interpreting CSV output")
    parser.add_argument(
        "--xwords", type=str, default="data/exclude_words.txt",
        help="Path to text file containing words to exclude. Defaults to "
             "data/exclude_words.txt")
    parser.add_argument(
        "--min_count", type=int, default=2,
        help="Minimum overall frequency of a word to be included in "
             "co-occurrence analysis.")
    parser.add_argument(
        "--min_records", type=int, default=3,
        help="Minimum number of records to include author (default 3).")
    args = parser.parse_args()

    xwords = load_xwords(args.xwords)
    files = glob(os.path.join(args.indir, "*.json"))
    debug = False
    if args.counts is None or args.matrix is None:
        debug = True

    # where all the action happens!
    all_words, by_author, by_document = process(
        files, xwords=xwords, debug=debug, min_records=args.min_records)

    # export indexed words
    intermed = [(count, word) for word, count in all_words.items()]
    intermed.sort(reverse=True)
    writer = csv.writer(args.index)
    writer.writerow(["word", "count", "index"])
    index = {}
    for idx, (count, word) in enumerate(intermed):
        if count < args.min_count:
            continue
        writer.writerow([word, count, idx])
        index.update({word: idx})

    if debug:
        sys.exit()  # print word index and quit

    # write counts by author
    json.dump(by_author, args.counts, indent=2)

    # export co-occurrences (context-term) as sparse matrix
    writer = csv.writer(args.matrix)
    for i, counts in enumerate(by_document):
        for word, count in counts.items():
            j = index.get(word, None)  # global indexing
            if j is not None:
                # document index, word index, count
                writer.writerow([i, j, count])
