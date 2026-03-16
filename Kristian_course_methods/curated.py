# -------------------------------------------------------
# FASTA filter that extracts real IDs from UniProt headers
# -------------------------------------------------------

# Your list of IDs to keep
KEEP = {
    id.strip()
    for id in """
YOUR ID1 FROM TREE
YOUR ID2 FROM TREE
...
""".strip().split()
}


def extract_id(header_line):
    """
    Extract real sequence ID from diverse FASTA headers.
    - UniProt format: >tr|A0A091MNE3|...
        → return "A0A091MNE3"
    - Simple format: >ITL03
        → return "ITL03"
    """
    header = header_line[1:].strip()

    if "|" in header:
        parts = header.split("|")
        return parts[1]  # accession field
    else:
        return header.split()[0]  # take up to first space


def filter_fasta(input_fasta, output_fasta):
    with open(input_fasta, "r") as infile, open(output_fasta, "w") as outfile:
        keep = False
        for line in infile:
            if line.startswith(">"):
                seq_id = extract_id(line)
                keep = seq_id in KEEP

                if keep:
                    outfile.write(line)

            else:
                if keep:
                    outfile.write(line)

    print("Done! Wrote all matching sequences to:", output_fasta)

filter_fasta("OWN1.faa", "curated.faa")