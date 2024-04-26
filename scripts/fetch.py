from Bio import Entrez
import argparse
import sys

description = """
Query PubMed for a given author name and retrieve title, abstract and keywords.
Write the results to a text file, one article per line.
"""


def extract(article):
	""" Extract abstract text, title and keywords from article record """
	# parse title and abstract text
	medline = article["MedlineCitation"]["Article"]
	title = medline["ArticleTitle"]
	abstract = medline.get("Abstract", None)
	abstext = ""
	if abstract:
		abstext = ' '.join([str(entry) for entry in abstract["AbstractText"]])
	
	# parse keywords
	keylist = article["MedlineCitation"]["KeywordList"]
	keywords = []
	for el in keylist:
		keywords.extend([str(kw) for kw in el])

	return f"{title} {abstext} {' '.join(keywords)}"


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
	args = parser.parse_args()
	
	Entrez.email = args.email
	articles = fetch(args.query, retmax=args.retmax)
	for article in articles:
		args.outfile.write(f"{extract(article)}\n")
