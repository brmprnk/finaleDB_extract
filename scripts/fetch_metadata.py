"""Fetch the FinaleDB metadata of a list of records and write it to a CSV file.

    python scripts/fetch_metadata.py IDS -o metadata.csv

IDS is a plain list of record IDs (EE87786, ...) or a CSV with an ID column (e.g. a sample sheet). One row
per requested ID is written, in the order of IDS. Records for which the FinaleDB API returns no result get
status=no_result and empty metadata: in the Jiang et al. 2015 study these are exactly the records whose hg38 fragment
file is empty (see README, "Records without hg38 data"). Optional step: the download needs only the IDs.

Columns: ID, status, sample_name, assay, readlen, instrument, seq_layout, fragNum (hg38), publication_name,
publication_doi, tissue, disease.
"""

import argparse
import csv
import sys
import time

import requests

API = "http://finaledb.research.cchmc.org/api/v1/seqrun"
COLUMNS = ["ID", "status", "sample_name", "assay", "readlen", "instrument", "seq_layout", "fragNum", "publication_name", "publication_doi", "tissue", "disease"]


def read_ids(path):
    """IDs from a plain list or from the ID column of a CSV; comment and blank lines skipped, duplicates dropped."""
    with open(path, newline="") as handle:
        lines = [line.strip() for line in handle if line.strip() and not line.lstrip().startswith("#")]
    if not lines:
        return []
    header = [h.strip().strip('"') for h in next(csv.reader([lines[0]]))]
    if "ID" in header:
        col = header.index("ID")
        ids = [row[col].strip() for row in csv.reader(lines[1:]) if row]
    else:
        ids = [line.split(",")[0].strip() for line in lines]
    return list(dict.fromkeys(i for i in ids if i))


def query(seq_id, session, retries=3):
    """The API result for one record, or None when the API returns no result."""
    for attempt in range(retries):
        try:
            response = session.get(API, params={"id": seq_id}, timeout=60)
            response.raise_for_status()
            results = response.json().get("results")
            if isinstance(results, list):
                return results[0] if results else None
            return results or None
        except (requests.RequestException, ValueError) as error:
            if attempt == retries - 1:
                raise RuntimeError(f"FinaleDB API request for {seq_id} failed: {error}") from error
            time.sleep(10 * (attempt + 1))
    return None


def to_row(seq_id, result):
    """One CSV row; the record's own ID is checked against the requested one."""
    if result is None or "sample" not in result:
        return {"ID": seq_id, "status": "no_result"}
    returned = f"EE{result.get('id')}"
    if returned != seq_id:
        raise RuntimeError(f"FinaleDB returned record {returned} for {seq_id}")
    seq_config = result.get("seqConfig", {})
    publication = result.get("publication", {})
    return {
        "ID": seq_id,
        "status": "ok",
        "sample_name": result["sample"].get("name"),
        "assay": result.get("assay"),
        "readlen": seq_config.get("readlen"),
        "instrument": seq_config.get("instrument"),
        "seq_layout": seq_config.get("seq_layout"),
        "fragNum": result.get("fragNum", {}).get("hg38"),
        "publication_name": publication.get("citeShort"),
        "publication_doi": publication.get("identifiers", {}).get("doi"),
        "tissue": result["sample"].get("tissue"),
        "disease": result["sample"].get("disease"),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("ids", help="plain ID list or CSV with an ID column")
    parser.add_argument("-o", "--output", required=True, help="output CSV")
    args = parser.parse_args()

    ids = read_ids(args.ids)
    print(f"Querying FinaleDB for {len(ids)} records")
    rows = []
    with requests.Session() as session:
        for seq_id in ids:
            rows.append(to_row(seq_id, query(seq_id, session)))
    with open(args.output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    missing = [row["ID"] for row in rows if row["status"] != "ok"]
    print(f"Wrote {len(rows)} rows to {args.output}; {len(missing)} without an API result" + (f": {' '.join(missing)}" if missing else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
