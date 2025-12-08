#!/usr/bin/env python3
"""
Build a metadata table for protein accessions using NCBI E-utilities.

Inputs:
  - One or more FASTA files (can be .gz) OR a text file with one accession per line.

Outputs:
  - TSV with columns:
    AccessionVersion, Id, Length, UpdateDate, TaxId, ScientificName, Title,
    Status, SeqType, Genome, Strain, Lineage

Usage:
  python entrez_metadata_from_accessions.py \
      --input virus_refseq_all_final_nr.fasta.gz \
      --out manifests/virus_refseq_manifest.tsv \
      --email elif.aksoy@bezmialem.edu.tr \
      --api-key $NCBI_API_KEY \
      --log run.log

Notes:
  - This *only* fetches metadata (ESummary + Taxonomy); no sequences are fetched.
  - Add --no-lineage to skip Taxonomy (faster).
  - Add --log <path> to write log messages to a file (still prints to stdout).
"""

import os, re, sys, gzip, time, csv, argparse
from typing import Iterable, List, Dict, Set
from Bio import Entrez

# ---------- Helpers ----------
_log_file = None

def log(msg: str):
    ts = time.strftime("[%Y-%m-%d %H:%M:%S]")
    log_line = f"{ts} {msg}\n"
    print(f"{ts} {msg}", flush=True)
    if _log_file is not None:
        _log_file.write(log_line)
        _log_file.flush()

def detect_is_text_list(path: str) -> bool:
    return os.path.splitext(path)[1].lower() in {".txt", ".list", ".acc", ".ids"}

# NEW: basic validation of protein accessions / IDs
_valid_refseq = re.compile(r"^[A-Z]{1,3}_\d+(?:\.\d+)?$")   # NP_123456.1, YP_, XP_, WP_ etc.
_valid_uniprot = re.compile(r"^[A-NR-Z0-9]{6,10}$")         # Uniprot-like
_valid_numeric = re.compile(r"^\d+$")                       # numeric UID

def is_valid_protein_id(s: str) -> bool:
    """
    Heuristic check for NCBI/UniProt-like protein identifiers.

    Rejects junk like 'AAA_28', '5AUM_D', '1Q3Z_A', etc.,
    which are typically PDB chain IDs or malformed headers.
    """
    s = s.strip()
    if not s:
        return False
    if _valid_numeric.match(s):
        return True
    if _valid_refseq.match(s):
        return True
    if _valid_uniprot.match(s):
        return True
    return False

def load_accessions_from_fasta(paths: List[str]) -> List[str]:
    accs: Set[str] = set()
    pat = re.compile(r"([A-Z]{1,3}_\d+(?:\.\d+)?)")  # e.g., NP_123456.1 / YP_ / WP_ / XP_
    for p in paths:
        opener = gzip.open if p.endswith(".gz") else open
        with opener(p, "rt") as fh:
            for line in fh:
                if line.startswith(">"):
                    # try regex first
                    m = pat.search(line)
                    token = None
                    if m:
                        token = m.group(1)
                    else:
                        # fallback: first token after ">"
                        token = line[1:].strip().split()[0]

                    # only keep if it's a valid protein id
                    if token and is_valid_protein_id(token):
                        accs.add(token)
                    else:
                        # optional: log once in a while for debugging
                        # (this won't spam too much on big datasets)
                        log(f"Skipping invalid/odd accession in FASTA header: {token!r}")
    return sorted(accs)

def load_accessions_from_text(path: str) -> List[str]:
    accs: Set[str] = set()
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if not s:
                continue
            if is_valid_protein_id(s):
                accs.add(s)
            else:
                log(f"Skipping invalid/odd accession from text list: {s!r}")
    return sorted(accs)

def chunker(seq: List[str], n: int) -> Iterable[List[str]]:
    for i in range(0, len(seq), n):
        yield seq[i:i+n]

def esummary_protein(ids: List[str], delay: float, retries: int) -> List[Dict]:
    """
    Fetch ESummary for protein IDs.

    Extra robustness:
      - If NCBI throws 'Invalid uid XXX', we drop XXX from the batch and retry
        instead of killing the whole run.
    """
    out: List[Dict] = []
    for batch in chunker(ids, 400):
        batch = list(batch)  # we may mutate it
        if not batch:
            continue

        attempt = 0
        while True:
            try:
                h = Entrez.esummary(db="protein", id=",".join(batch), retmode="xml")
                recs = Entrez.read(h)
                h.close()
                out.extend(recs)
                time.sleep(delay)
                break
            except Exception as e:
                attempt += 1
                msg = str(e)

                # Try to detect and remove specific invalid uid from the batch
                m = re.search(r"Invalid uid\s+(\S+)", msg)
                if m:
                    bad_uid = m.group(1)
                    if bad_uid in batch:
                        log(f"Entrez ESummary reported invalid uid {bad_uid!r}; dropping from batch and retrying.")
                        batch = [x for x in batch if x != bad_uid]
                        if not batch:
                            log("All IDs in this batch were removed as invalid; skipping batch.")
                            break
                        # do *not* count this as a real retry; continue immediately
                        continue

                if attempt > retries:
                    log(f"ESummary failed after {retries} retries; last error: {msg}")
                    raise
                log(f"ESummary retry {attempt}/{retries} after error: {msg}")
                time.sleep(3 * attempt)
    return out

def taxonomy_lineage(taxids: List[str], delay: float, retries: int) -> Dict[str, Dict[str, str]]:
    lineage_map: Dict[str, Dict[str, str]] = {}
    for batch in chunker(taxids, 400):
        attempt = 0
        while True:
            try:
                h = Entrez.efetch(db="taxonomy", id=",".join(batch), retmode="xml")
                taxa = Entrez.read(h)
                h.close()
                for t in taxa:
                    tid = str(t.get("TaxId", ""))
                    lineage_map[tid] = {
                        "ScientificName": t.get("ScientificName", ""),
                        "Lineage": t.get("Lineage", ""),
                    }
                time.sleep(delay)
                break
            except Exception as e:
                attempt += 1
                if attempt > retries:
                    log(f"Taxonomy failed after {retries} retries; last error: {e}")
                    raise
                log(f"Taxonomy retry {attempt}/{retries} after error: {e}")
                time.sleep(3 * attempt)
    return lineage_map

# ---------- Main ----------
def main():
    global _log_file
    
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, nargs="+",
                    help="FASTA(.gz) and/or a text file of accessions")
    ap.add_argument("--out", required=True, help="Output TSV path")
    ap.add_argument("--email", required=True, help="Your email for NCBI E-utilities")
    ap.add_argument("--api-key", default=None, help="NCBI API key (optional)")
    ap.add_argument("--no-lineage", action="store_true", help="Skip taxonomy lineage fetch")
    ap.add_argument("--delay", type=float, default=None,
                    help="Seconds between calls (auto-sets by API key if not provided)")
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--log", default=None, help="Path to log file (optional)")
    args = ap.parse_args()

    # Open log file if provided
    if args.log:
        os.makedirs(os.path.dirname(args.log) or ".", exist_ok=True)
        _log_file = open(args.log, "w")
        log(f"Logging to: {args.log}")

    try:
        Entrez.email = args.email
        if args.api_key:
            Entrez.api_key = args.api_key

        delay = args.delay
        if delay is None:
            delay = 0.10 if args.api_key else 0.34  # polite defaults

        # Collect accessions
        fasta_paths, text_paths = [], []
        for p in args.input:
            if detect_is_text_list(p):
                text_paths.append(p)
            else:
                fasta_paths.append(p)

        accs: List[str] = []
        if fasta_paths:
            log(f"Scanning FASTA headers from: {', '.join(fasta_paths)}")
            accs.extend(load_accessions_from_fasta(fasta_paths))
        for t in text_paths:
            log(f"Reading accessions from: {t}")
            accs.extend(load_accessions_from_text(t))

        # Deduplicate
        accs = sorted(set(accs))
        if not accs:
            log("No (valid) accessions found. Exiting.")
            sys.exit(1)

        log(f"Total unique *valid* accessions before ESummary: {len(accs)}")

        # ESummary
        log("Fetching ESummary metadata (protein)…")
        summaries = esummary_protein(accs, delay=delay, retries=args.retries)

        # Build rows + collect TaxIds
        rows = []
        taxids = set()
        for d in summaries:
            # ESummary keys are somewhat variable across records; be defensive
            av = d.get("AccessionVersion", d.get("Caption", ""))
            row = {
                "AccessionVersion": av,
                "Id": str(d.get("Id", "")),
                "Length": str(d.get("Length", "")),
                "UpdateDate": d.get("UpdateDate", ""),
                "TaxId": str(d.get("TaxId", "")),
                "ScientificName": d.get("ScientificName", ""),  # sometimes present
                "Title": d.get("Title", ""),
                "Status": d.get("Status", ""),
                "SeqType": d.get("MoleculeType", d.get("MolType", "")),
                "Genome": d.get("Genome", ""),
                "Strain": d.get("SubName", ""),   # may be empty
                "Lineage": "",                    # fill after taxonomy
            }
            rows.append(row)
            if row["TaxId"]:
                taxids.add(row["TaxId"])

        # Taxonomy lineage (optional)
        lineage_map = {}
        if not args.no_lineage and taxids:
            log(f"Fetching Taxonomy lineage for {len(taxids)} TaxIds…")
            lineage_map = taxonomy_lineage(sorted(taxids), delay=delay, retries=args.retries)

        for r in rows:
            tid = r["TaxId"]
            if tid and tid in lineage_map:
                # If ESummary lacked ScientificName, fill from taxonomy
                if not r["ScientificName"]:
                    r["ScientificName"] = lineage_map[tid].get("ScientificName", "")
                r["Lineage"] = lineage_map[tid].get("Lineage", "")

        # Write TSV
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        fields = ["AccessionVersion","Id","Length","UpdateDate","TaxId",
                  "ScientificName","Title","Status","SeqType","Genome","Strain","Lineage"]
        with open(args.out, "w", newline="") as fh:
            w = csv.DictWriter(fh, delimiter="\t", fieldnames=fields)
            w.writeheader()
            for r in rows:
                w.writerow({k: r.get(k, "") for k in fields})

        log(f"Done. Wrote {len(rows)} records → {args.out}")
    finally:
        if _log_file is not None:
            _log_file.close()

if __name__ == "__main__":
    main()
