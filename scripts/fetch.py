from Bio import Entrez
import argparse
import sys
import json

# TODO: convert .bib file to JSON?
description = """
Query PubMed for a given author name and retrieve title, abstract and keywords.
Write the results to a JSON file.
"""


def extract(article, max_words=1000):
    """
    Extract abstract text, title and keywords from article record
    :param article:
    :param max_words:  int, if abstract word count exceeds this number, discard
                       the record
    """
    # parse title and abstract text
    medline = article["MedlineCitation"]
    article = medline["Article"]
    
    # try to retrieve year of publication
    try:
        year = article["Journal"]['JournalIssue']['PubDate']['Year']
    except KeyError:
        # backup
        try:
            year = article["ArticleDate"][0]['Year']
        except:
            print(article)
            raise
        
    authors = []
    for aut in article["AuthorList"]:
        if "CollectiveName" in aut:
            authors.append(aut["CollectiveName"])
        else:
            names = [aut[key] for key in ["ForeName", "LastName"] if key in aut]
            authors.append(' '.join(names))
    title = article["ArticleTitle"]
    if title.startswith("Author Correction") or title.startswith("Correction:"):
        return None
    abstract = article.get("Abstract", None)
    abstext = ""
    if abstract:
        abstext = ' '.join([str(entry) for entry in abstract["AbstractText"]])
        if len(abstext.split()) > max_words:
            sys.stderr.write(f"Error: excessive word count in \"{title}\", skipping\n")
            return None

    # parse keywords
    keylist = medline["KeywordList"]
    keywords = []
    for el in keylist:
        keywords.extend([str(kw) for kw in el])

    return {'pmid': str(medline["PMID"]), 'year': year, 'authors': authors, 'title': title,
            'abstract': abstext, 'keywords': keywords}


def fetch(query, retmax=100):
    """ Handles PubMed API transaction.  query should be author name """
    response = Entrez.esearch(db="pubmed", term=query, retmax=retmax)
    result = Entrez.read(response)
    hits = Entrez.efetch(db="pubmed", rettype="abstract", id=result["IdList"])
    records = Entrez.read(hits)
    articles = records["PubmedArticle"]
    return articles


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("query", type=str, help="Author query")
    parser.add_argument("--email", type=str, default="apoon42@uwo.ca",
                        help="e-mail address for Entrez transactions")
    parser.add_argument("-o", "--outfile", type=argparse.FileType('w'),
                     default=sys.stdout, help="Output file, defaults to stdout.")
    parser.add_argument("--retmax", type=int, default=20, 
                     help="Maximum number of records to retrieve (default 20).")
    parser.add_argument("--max_words", type=int, default=1000,
                        help="Limit abstract word count, to avoid problematic recods (default 1000).")
    args = parser.parse_args()
    
    Entrez.email = args.email
    articles = fetch(args.query, retmax=args.retmax)

    records = [extract(art, max_words=args.max_words) for art in articles]
    records = [r for r in records if r is not None]
    sys.stderr.write(f"Retrieved {len(records)} records (limit {args.retmax})\n")
    
    json.dump(records, args.outfile, indent=2)

