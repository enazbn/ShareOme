#!/usr/bin/env python3
"""
Clean protein FASTA files before kmerslicing, with robust Entrez validation.

Goals:
  - Remove only *truly* invalid / non-protein accessions (AAA_28, BBA_00009, nuccore-only UIDs, etc.).
  - NEVER drop an entire batch because of one bad ID.
  - Handle 'Otherdb uid="47118297" db="nuccore"' cleanly.
  - Avoid infinite retry loops.

Strategy:
  1) Pattern-based filtering (optional / minimal): remove obvious PDB chains like 1ABC_A.
  2) Entrez-based validation:
       - First try batches (e.g. 400 IDs).
       - If batch fails with 'Invalid uid X' and X is in the batch → mark X invalid, remove, retry batch.
       - For any other persistent batch error (e.g. Otherdb weirdness) → fall back to per-ID validation for that batch.
         No batch is ever discarded wholesale.

Output:
  - A cleaned FASTA (.fasta or .fasta.gz) containing only accessions that:
      * pass minimal pattern checks, and
      * are accepted by NCBI ESummary as protein records.

Usage:
  python clean_fasta_for_shareome.py \
      --input virus_refseq_all_final_nr \
      --out virus_refseq_all_final_nr.cleaned.fasta \
      --email YOUR_EMAIL \
      --api-key $NCBI_API_KEY \
      --log logs/clean_fasta.log
"""

import os
import re
import sys
import gzip
import time
import argparse
from typing import List, Set
from Bio import Entrez
import socket
socket.setdefaulttimeout(60)  # or 60 if you prefer

# ---------- Logging ----------
_log_file = None

def log(msg: str):
    ts = time.strftime("[%Y-%m-%d %H:%M:%S]")
    line = f"{ts} {msg}"
    print(line, flush=True)
    if _log_file is not None:
        _log_file.write(line + "\n")
        _log_file.flush()

# ---------- ID handling ----------

def is_valid_protein_id_local(s: str) -> bool:
    """
    Minimal local validation: only reject obvious junk like PDB chains.

    Pattern rejected:
      digit + 3 alnum + '_' + 1 alnum  => classic PDB chain IDs (1Q3Z_A).

    Everything else is allowed here. Entrez will decide what's truly invalid.
    """
    s = s.strip()
    if not s:
        return False

    # PDB chain-like ID, e.g. 1Q3Z_A
    if re.match(r"^[0-9][A-Za-z0-9]{3}_[A-Za-z0-9]$", s):
        return False

    return True

def extract_accession_from_header(line: str) -> str:
    """
    Extract the accession from a FASTA header line.

    Logic:
      - Try to grab RefSeq-style ID like NP_123456.1 (NP_/XP_/YP_/WP_ etc.).
      - Otherwise, use the first token after '>'.
    """
    line = line.rstrip("\n")
    if not line.startswith(">"):
        return ""

    pat = re.compile(r"([A-Z]{1,3}_\d+(?:\.\d+)?)")
    m = pat.search(line)
    if m:
        return m.group(1)

    # Fallback: first token after '>'
    token = line[1:].strip().split()[0]
    return token

# ---------- Entrez helpers ----------

def chunker(seq: List[str], n: int):
    for i in range(0, len(seq), n):
        yield seq[i:i+n]

def esummary_batch_ok(batch: List[str], delay: float) -> bool:
    """
    Helper to call ESummary once for a batch and return True if it works.

    Used only in the per-ID fallback.
    """
    h = Entrez.esummary(db="protein", id=",".join(batch), retmode="xml")
    _ = Entrez.read(h)
    h.close()
    time.sleep(delay)
    return True

def validate_ids_with_entrez(
    ids: List[str],
    email: str,
    api_key: str = None,
    delay: float = None,
    retries: int = 3,
    batch_size: int = 400,
) -> Set[str]:
    """
    Entrez-based validation of accessions.

    Returns:
      invalid_ids: set of IDs that failed validation.

    Behavior:
      - Try batch ESummary (size=batch_size).
      - If ESummary OK → all IDs in batch are considered valid.
      - If error with "Invalid uid X" and X in batch:
            mark X invalid, drop from batch, retry the batch.
      - For any other persistent batch error (e.g. Otherdb weirdness):
            fall back to per-ID validation for that batch:
              * For each accession:
                    try ESummary(id=that_one)
                    if it raises → mark that accession invalid.
            → no batch is ever discarded wholesale.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    if delay is None:
        delay = 0.10 if api_key else 0.34  # polite defaults

    invalid_ids: Set[str] = set()
    total = len(ids)
    processed = 0

    log(f"Entrez validation of {total} IDs (batch_size={batch_size})...")

    for original_batch in chunker(ids, batch_size):
        batch = list(original_batch)
        if not batch:
            continue

        attempt = 0
        while True:
            try:
                h = Entrez.esummary(db="protein", id=",".join(batch), retmode="xml")
                _ = Entrez.read(h)
                h.close()
                processed += len(batch)
                log(f"  ESummary OK for batch ({processed}/{total})")
                time.sleep(delay)
                break  # batch done, go to next batch

            except Exception as e:
                msg = str(e)

                # Case 1: "Invalid uid XXX" we can map directly
                m_invalid = re.search(r"Invalid uid\s+(\S+)", msg)
                if m_invalid:
                    bad_uid = m_invalid.group(1)
                    if bad_uid in batch:
                        invalid_ids.add(bad_uid)
                        log(f"  Entrez: Invalid uid {bad_uid!r}; dropping from batch and retrying.")
                        batch = [x for x in batch if x != bad_uid]
                        if not batch:
                            log("  All IDs in this batch were invalid uids; skipping batch.")
                            break
                        # handled, do NOT count as retry
                        continue

                # For anything else (e.g. Otherdb uid="47118297") we don't trust batch-level info.
                attempt += 1
                if attempt > retries:
                    # Fall back to per-ID validation for this batch
                    log(f"  Batch-level ESummary failed after {retries} retries with error: {msg}")
                    log("  Falling back to per-ID validation for this batch.")
                    for acc in batch:
                        try:
                            h2 = Entrez.esummary(db="protein", id=acc, retmode="xml")
                            _ = Entrez.read(h2)
                            h2.close()
                            time.sleep(delay)
                        except Exception as e2:
                            invalid_ids.add(acc)
                            log(f"    Per-ID: {acc!r} INVALID ({e2})")
                    processed += len(batch)
                    log(f"  Per-ID fallback complete for batch. Progress: ({processed}/{total})")
                    break  # move to next outer batch

                log(f"  ESummary retry {attempt}/{retries} after error: {msg}")
                time.sleep(3 * attempt)

    log(f"Entrez validation complete. Invalid IDs found: {len(invalid_ids)}")
    return invalid_ids

# ---------- FASTA cleaning ----------

def collect_accessions_from_fasta(path: str) -> List[str]:
    """
    Scan a FASTA(.gz) file and collect all accession IDs from headers.
    """
    opener = gzip.open if path.endswith(".gz") else open
    accs: Set[str] = set()
    headers = 0
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                headers += 1
                acc = extract_accession_from_header(line)
                if acc:
                    accs.add(acc)
    log(f"Scanned {headers} headers in {path}; found {len(accs)} unique raw accessions.")
    return sorted(accs)

def write_cleaned_fasta(
    in_path: str,
    out_path: str,
    bad_ids: Set[str]
):
    """
    Rewrite FASTA, skipping any record whose accession is in bad_ids.
    """
    in_open = gzip.open if in_path.endswith(".gz") else open
    out_open = gzip.open if out_path.endswith(".gz") else open

    kept = 0
    skipped = 0

    with in_open(in_path, "rt") as fin, out_open(out_path, "wt") as fout:
        write_current = True

        for line in fin:
            if line.startswith(">"):
                acc = extract_accession_from_header(line)
                if acc in bad_ids:
                    write_current = False
                    skipped += 1
                else:
                    write_current = True
                    kept += 1
                    fout.write(line)
            else:
                if write_current:
                    fout.write(line)

    log("FASTA cleaning complete.")
    log(f"  Records kept   : {kept}")
    log(f"  Records removed: {skipped}")
    log(f"Cleaned FASTA written to: {out_path}")

# ---------- Main ----------

def main():
    global _log_file

    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True,
                    help="Input protein FASTA (.fasta or .fasta.gz)")
    ap.add_argument("--out", required=True,
                    help="Output cleaned FASTA (.fasta or .fasta.gz)")
    ap.add_argument("--email", required=True,
                    help="Your email for NCBI E-utilities")
    ap.add_argument("--api-key", default=None,
                    help="NCBI API key (optional, but recommended)")
    ap.add_argument("--delay", type=float, default=None,
                    help="Seconds between Entrez calls (auto-set if not provided)")
    ap.add_argument("--retries", type=int, default=3,
                    help="Max retries for batch-level Entrez calls before per-ID fallback")
    ap.add_argument("--log", default=None,
                    help="Log file path (optional)")
    args = ap.parse_args()

    # Set up logging
    if args.log:
        os.makedirs(os.path.dirname(args.log) or ".", exist_ok=True)
        _log_file = open(args.log, "w")
        log(f"Logging to: {args.log}")

    try:
        if not os.path.exists(args.input):
            log(f"ERROR: Input FASTA does not exist: {args.input}")
            sys.exit(1)

        log(f"Input FASTA : {args.input}")
        log(f"Output FASTA: {args.out}")

        # 1) Collect all raw accessions from FASTA
        all_accs = collect_accessions_from_fasta(args.input)

        # 2) Minimal pattern-based invalid IDs (PDB-like, etc.)
        bad_pattern_ids = {a for a in all_accs if not is_valid_protein_id_local(a)}
        good_for_entrez = [a for a in all_accs if a not in bad_pattern_ids]

        log(f"Pattern-based invalid IDs (e.g. PDB-like): {len(bad_pattern_ids)}")
        if bad_pattern_ids:
            log(f"Example bad pattern IDs: {list(sorted(bad_pattern_ids))[:10]}")

        # 3) Entrez validation for remaining IDs
        bad_entrez_ids = validate_ids_with_entrez(
            good_for_entrez,
            email=args.email,
            api_key=args.api_key,
            delay=args.delay,
            retries=args.retries
        )

        # 4) Combine all bad IDs
        bad_ids = set(bad_pattern_ids) | set(bad_entrez_ids)
        log(f"Total bad IDs to be removed from FASTA: {len(bad_ids)}")

        # 5) Write cleaned FASTA
        os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
        write_cleaned_fasta(args.input, args.out, bad_ids)

    finally:
        if _log_file is not None:
            _log_file.close()

if __name__ == "__main__":
    main()
