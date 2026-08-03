import os
import pandas as pd
from Bio import SeqIO

base_dir = "/media/ramesh/LGI/myb-01072026/NLR–NAC–MYB-Axis"
out_dir = os.path.join(base_dir, "analysis_results")
os.makedirs(out_dir, exist_ok=True)

# 1. Load Name Mappings
def load_mapping(file_path, family_name):
    mappings = []
    if not os.path.exists(file_path):
        print(f"Warning: {file_path} not found")
        return mappings
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 2:
                orig_id, new_id = parts[0], parts[1]
                # Normalize ID (e.g., AH012382-RA or AH012382)
                gene_id = orig_id.split('.')[0]
                mappings.append({
                    'Original_ID': orig_id,
                    'Gene_ID': gene_id,
                    'New_ID': new_id,
                    'Family': family_name
                })
    return mappings

myb_map = load_mapping(os.path.join(base_dir, "myb-old-new-names.txt"), "MYB")
nac_map = load_mapping(os.path.join(base_dir, "nac-old-new-names.txt"), "NAC")
nlr_map = load_mapping(os.path.join(base_dir, "nlr-old-new-names.txt"), "NLR")

df_meta = pd.DataFrame(myb_map + nac_map + nlr_map)
meta_path = os.path.join(out_dir, "master_target_genes_meta.csv")
df_meta.to_csv(meta_path, index=False)

print(f"Loaded total metadata entries: {len(df_meta)}")
print(df_meta['Family'].value_counts())

# 2. Extract Protein Sequences for all 311 genes
prot_out = os.path.join(out_dir, "target_proteins_combined.fasta")
fasta_files = [
    (os.path.join(base_dir, "myb-sorted-renamed.fasta"), "MYB"),
    (os.path.join(base_dir, "nac-sorted-renamed.fasta"), "NAC"),
    (os.path.join(base_dir, "nlr-sorted-renamed.fasta"), "NLR")
]

records = []
for fpath, fam in fasta_files:
    if os.path.exists(fpath):
        for rec in SeqIO.parse(fpath, "fasta"):
            records.append(rec)

SeqIO.write(records, prot_out, "fasta")
print(f"Wrote {len(records)} protein sequences to {prot_out}")

# 3. Extract 2.0 kb Promoter Sequences
prom_file = os.path.join(base_dir, "Ah-Promoter.fasta")
prom_out = os.path.join(out_dir, "target_promoters_2kb.fasta")

target_orig_ids = set(df_meta['Original_ID']).union(set(df_meta['Gene_ID']))
prom_records = []

if os.path.exists(prom_file):
    for rec in SeqIO.parse(prom_file, "fasta"):
        # Match header ID
        rec_id = rec.id.split()[0]
        rec_base = rec_id.split('.')[0]
        if rec_id in target_orig_ids or rec_base in target_orig_ids:
            prom_records.append(rec)

SeqIO.write(prom_records, prom_out, "fasta")
print(f"Extracted {len(prom_records)} promoter sequences to {prom_out}")
