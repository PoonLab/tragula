import argparse
import json
import nltk
import os
from glob import glob
import sys

description = """
Tokenize article text for all author-specific JSON files in a given 
folder.  Export word co-occurrence as sparse matrix, and author 
word frequency vectors as a second file.
https://www.ibm.com/docs/en/wca/3.5.0?topic=analytics-part-speech-tag-sets
"""

xtags = {"IN", "CD", "DT", "CC", "TO", ":", "(", ")", ",", "."}
xwords = {
    "we", "is", "are", "were", "that", "have", "can", "which", "may",
    "<", ">", "=", "%", "was", "using", "be", "more", "not", "our", "had",
    "however", "p", "use", "their", "i", "used", "it", "has", "such", "t",
    "other", "been", "n", "''", "``", "including", "its", "most", "when", "+",
    "suggest", "significantly", "significant", "different", "show", "who"
}


def process(files,
            ):
    """
    Enumerate word co-occurrence across all texts, and track
    author-specific word frequencies.  Remove the following from English tag set:
    IN = preposition / subordinating conjunction
    CD = cardinal number
    DT = determiner
    CC = coordinating conjunction\
    TO = preposition "to"
    """
    by_author = {}
    corpus = {}
    for f in files:
        author = os.path.basename(f).split('.')[0]
        if author in corpus:
            sys.stderr.write(f"ERROR: duplicate author {author} detected!")
            sys.exit()
        by_author.update({author: {}})

        with open(f) as handle:
            records = json.load(handle)
        for r in records:
            text = f"{r['title']} {r['abstract']} {' '.join(r['keywords'])}"
            tokens = nltk.word_tokenize(text)
            tagged = nltk.pos_tag(tokens)  # part-of-speech tagger
            filtered = [word.lower() for word, tag in tagged
                        if tag not in xtags and word.lower() not in xwords]
            for word in filtered:
                if word not in corpus:
                    corpus.update({word: 0})
                corpus[word] += 1
                if word not in by_author[author]:
                    by_author[author].update({word: 0})
                by_author[author][word] += 1

    intermed = [(count, word) for word, count in corpus.items()]
    intermed.sort(reverse=True)
    print(intermed[:100])


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description)
    parser.add_argument("indir", type=str, help="Path to folder containing JSON files.")
    args = parser.parse_args()

    files = glob(os.path.join(args.indir, "*.json"))
    process(files)
