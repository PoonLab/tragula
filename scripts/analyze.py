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
word frequency vectors as a second file.  If no output files are 
specified, then the script will write all word frequencies to stdout 
in CSV format.
"""

# tags to exclude
# https://www.ibm.com/docs/en/wca/3.5.0?topic=analytics-part-speech-tag-sets
xtags = {
    "IN",  # preposition / subordinating conjunction
    "CD",  # cardinal number
    "DT",  # determiner
    "CC",  # coordinating conjunction
    "TO",  # preposition "to"
    ":", "(", ")", ",", "."  # punctuation
}

# words to exclude
xwords = {
    "", "we", "is", "are", "were", "that", "have", "can", "which", "may", "also",
    "a", "<", ">", "=", "%", "was", "using", "be", "more", "not", "our", "had",
    "however", "p", "use", "their", "i", "/i", "used", "it", "has", "such", "t",
    "other", "been", "n", "''", "``", "including", "its", "most", "when", "+",
    "suggest", "significantly", "significant", "different", "show", "who",
    "always", "study", "increased", "change", "result", "associated", "method",
    "number", "high", "increase", "compared", "response", "showed",
    "identified", "further", "found", "well", "examined", "specific", "novel",
    "following", "sup", "/sup", "observed", "important", "there", "research",
    "here", "evidence", "year", "month", "furthermore", "demonstrated",
    "reduced", "several", "multiple", "measured", "recent", "new", "known",
    "only", "[", "]", "review", "common", "revealed", "large", "assessed",
    "approach", "shown", "indicate", "based", "no", "'", "decreased",
    "caused", "due", "determine", "could", "previously", "often", "they",
    "will", "same", "studied", "previous", "many", "understanding",
    "lower", "included", "improved", "underlying", "then", "evaluate",
    "required", "remains", "wa"
}

wordnet = WordNetLemmatizer()  # convert plural to singular
numeric = re.compile("^-?[0-9]+\.?[0-9]*$")


def process(files, debug=False):
    """
    Enumerate word co-occurrence across all texts, and track
    author-specific word frequencies.
    :param files:  list, paths to JSON files containing
    :param debug:  bool, get global word counts only
    :returns:
      - by_author - dict, word counts by author
      - cooccur - dict, stores numbers of texts in which both words appear
    """
    by_author = {}
    cooccur = {}
    all_words = {}
    for f in files:
        author = os.path.basename(f).split('.')[0]
        if author in by_author:
            sys.stderr.write(f"ERROR: duplicate author {author} detected!")
            sys.exit()
        by_author.update({author: {}})

        with open(f) as handle:
            records = json.load(handle)

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

            # increment word counts by author
            if debug:
                for word in filtered:
                    if word not in all_words:
                        all_words.update({word: 0})
                    all_words[word] += 1
                continue  # skip next step
            else:
                for word in filtered:
                    if word not in by_author[author]:
                        by_author[author].update({word: 0})
                    by_author[author][word] += 1

            # record co-occurrences
            intermed = list(set(filtered))
            intermed.sort()  # w1 < w2
            for i in range(len(intermed)-1):
                w1 = intermed[i]
                if w1 not in cooccur:
                    cooccur.update({w1: {}})
                for j in range(i+1, len(intermed)):
                    w2 = intermed[j]
                    if w2 not in cooccur[w1]:
                        cooccur[w1].update({w2: 0})
                    cooccur[w1][w2] += 1

    if debug:
        return all_words, None
    else:
        return by_author, cooccur


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description)
    parser.add_argument("indir", type=str, help="Path to folder containing JSON files.")
    parser.add_argument("--counts", type=argparse.FileType('w'), required=False,
                        help="Path to write by-author word counts (JSON)")
    parser.add_argument("--matrix", type=argparse.FileType('w'), required=False,
                        help="Path to write co-occurrence sparse matrix (CSV)")
    parser.add_argument("--index", type=argparse.FileType('w'), required=False,
                        help="Path to write indexed words (for interpreting CSV)")
    args = parser.parse_args()

    files = glob(os.path.join(args.indir, "*.json"))

    if args.counts is None or args.matrix is None:
        counts, _ = process(files, debug=True)
        intermed = [(count, word) for word, count in counts.items()]
        intermed.sort(reverse=True)
        for count, word in intermed:
            print(f"{word},{count}")
    else:
        by_author, cooccur = process(files, debug=False)
        json.dump(by_author, args.counts, indent=2)

        # collect all words
        lex = set()
        for w1, rows in cooccur.items():
            lex.add(w1)
            for w2, _ in rows.items():
                lex.add(w2)
        lex = sorted(lex)
        index = dict([(w, i) for i, w in enumerate(lex)])
        for w in lex:
            args.index.write(f"{w}\n")

        # export co-occurrences as sparse matrix
        writer = csv.writer(args.matrix)
        for w1, rows in cooccur.items():
            i = index[w1]
            for w2, count in rows.items():
                if count == 1:
                    continue
                j = index[w2]
                writer.writerow([i, j, count])
