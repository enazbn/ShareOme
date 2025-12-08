from Bio import Entrez
import time
import os
import gzip
import sys

# --- USER SETTINGS ---
Entrez.email = "elif.aksoy@bezmialem.edu.tr"  # CHANGE THIS to your own email
SAVE_DIR = "human_refseq_proteins"  # Directory to save downloaded files
BATCH_SIZE = 500  # Number of records to fetch per batch
RETRY_LIMIT = 3  # Number of retries on failure
DELAY = 0.5  # Delay between requests (in seconds)

# --- NCBI QUERY ---
QUERY = "txid9606[Organism]"

# --- Helper function for timestamped logging ---
def timestamp():
    return time.strftime("[%Y-%m-%d %H:%M:%S]")

# --- Ensure save directory exists ---
os.makedirs(SAVE_DIR, exist_ok=True)

# --- Start search ---
print(f"{timestamp()} Starting Human RefSeq Protein Download")
print(f"{timestamp()} Query used: {QUERY}")
handle = Entrez.esearch(db="protein", term=QUERY, usehistory="y")
search_results = Entrez.read(handle)
handle.close()

count = int(search_results["Count"])
webenv = search_results["WebEnv"]
query_key = search_results["QueryKey"]

print(f"{timestamp()} Found {count} protein records. Starting download...")

# --- Fetch sequences in batches ---
for start in range(0, count, BATCH_SIZE):
    end = min(count, start + BATCH_SIZE)
    print(f"{timestamp()} Fetching records {start + 1} to {end}...")

    success = False
    retries = 0
    while not success and retries < RETRY_LIMIT:
        try:
            fetch_handle = Entrez.efetch(
                db="protein",
                rettype="fasta",
                retmode="text",
                retstart=start,
                retmax=BATCH_SIZE,
                webenv=webenv,
                query_key=query_key,
            )
            data = fetch_handle.read()
            fetch_handle.close()

            filename = os.path.join(SAVE_DIR, f"human_refseq_{start+1}_{end}.fasta.gz")
            with gzip.open(filename, "wt") as out_handle:
                out_handle.write(data)

            success = True
        except Exception as e:
            print(f"{timestamp()} Error fetching batch {start}-{end}: {e}")
            retries += 1
            time.sleep(5)

    if not success:
        print(f"{timestamp()} Failed to fetch batch {start}-{end} after {RETRY_LIMIT} retries.")
        sys.exit(1)

    time.sleep(DELAY)

print(f"{timestamp()} Download complete!")
