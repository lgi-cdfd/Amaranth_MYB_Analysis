import os
import re
import pandas as pd
from Bio import SeqIO

base_dir = "/media/ramesh/LGI/myb-01072026/NLR–NAC–MYB-Axis"
out_dir = os.path.join(base_dir, "analysis_results")
os.makedirs(out_dir, exist_ok=True)

# 1. Load Promoter Sequences
prom_file = os.path.join(out_dir, "target_promoters_2kb.fasta")
records = list(SeqIO.parse(prom_file, "fasta"))

meta_path = os.path.join(out_dir, "master_target_genes_meta.csv")
df_meta = pd.read_csv(meta_path)

# Map header to New_ID and Family
id_map = {}
for idx, row in df_meta.iterrows():
    id_map[row['Original_ID']] = (row['New_ID'], row['Family'])
    id_map[row['Gene_ID']] = (row['New_ID'], row['Family'])

# Defined Cis-Regulatory Motifs
motifs = {
    'NAC_binding_site': [r'CGT[GA]', r'CACG', r'CATGT'],
    'MYB_binding_site': [r'[AT]AACCA', r'AACGG', r'CAGTTA', r'CCGTT'],
    'W_box_Immune': [r'TTGAC'],
    'ABRE_Drought': [r'ACGTG']
}

motif_counts = []
edges = []

for rec in records:
    rec_id = rec.id.split()[0]
    rec_base = rec_id.split('.')[0]
    
    new_id, fam = id_map.get(rec_id, id_map.get(rec_base, (rec_id, "Unknown")))
    seq_str = str(rec.seq).upper()
    
    row = {'Gene_ID': new_id, 'Family': fam, 'Original_ID': rec_id}
    
    for mname, mpatterns in motifs.items():
        total_hits = 0
        for pat in mpatterns:
            matches = re.findall(pat, seq_str)
            total_hits += len(matches)
        row[mname] = total_hits
        
        # Build regulatory edges (Target = this gene's promoter)
        if mname == 'NAC_binding_site' and total_hits >= 2:
            edges.append({'Source_Type': 'NAC_TF_Family', 'Target_Gene': new_id, 'Target_Family': fam, 'Interaction': 'Promoter_NAC_Motif', 'Weight': total_hits})
        elif mname == 'MYB_binding_site' and total_hits >= 2:
            edges.append({'Source_Type': 'MYB_TF_Family', 'Target_Gene': new_id, 'Target_Family': fam, 'Interaction': 'Promoter_MYB_Motif', 'Weight': total_hits})
            
    motif_counts.append(row)

df_motifs = pd.DataFrame(motif_counts)
motif_out = os.path.join(out_dir, "promoter_cis_regulatory_motifs_2kb.csv")
df_motifs.to_csv(motif_out, index=False)

df_edges = pd.DataFrame(edges)
edge_out = os.path.join(out_dir, "cis_regulatory_network_edges.csv")
df_edges.to_csv(edge_out, index=False)

print(f"Scanned {len(motif_counts)} promoters for key cis-elements.")
print("\n--- Average Motif Density per Family ---")
print(df_motifs.groupby('Family')[['NAC_binding_site', 'MYB_binding_site', 'W_box_Immune', 'ABRE_Drought']].mean().round(2))

# 2. Build Protein-Protein Interaction (PPI) Co-clustering Edges
# Load WGCNA modules
modules_path = os.path.join(out_dir, "wgcna_coexpression_modules_310_genes.csv")
df_mod = pd.read_csv(modules_path)

ppi_edges = []
# Connect high-kME hub NLRs to NACs and MYBs in the same co-expression module
for mod_id in df_mod['Module_ID'].unique():
    mod_sub = df_mod[df_mod['Module_ID'] == mod_id]
    
    nlrs = mod_sub[mod_sub['Family'] == 'NLR']['Gene_ID'].tolist()
    nacs = mod_sub[mod_sub['Family'] == 'NAC']['Gene_ID'].tolist()
    mybs = mod_sub[mod_sub['Family'] == 'MYB']['Gene_ID'].tolist()
    
    # Create NLR -> NAC and NAC -> MYB cross-family edges
    for nlr in nlrs[:5]: # Top 5 NLR hubs
        for nac in nacs[:5]: # Top 5 NAC hubs
            ppi_edges.append({'Source': nlr, 'Target': nac, 'Source_Fam': 'NLR', 'Target_Fam': 'NAC', 'Module': mod_id, 'Edge_Type': 'NLR-NAC_CoResponse'})
            
    for nac in nacs[:5]:
        for myb in mybs[:5]:
            ppi_edges.append({'Source': nac, 'Target': myb, 'Source_Fam': 'NAC', 'Target_Fam': 'MYB', 'Module': mod_id, 'Edge_Type': 'NAC-MYB_Regulation'})

df_ppi = pd.DataFrame(ppi_edges)
ppi_out = os.path.join(out_dir, "nlr_nac_myb_interactome_edges.csv")
df_ppi.to_csv(ppi_out, index=False)

print(f"\nGenerated {len(df_ppi)} interactome edges across the NLR-NAC-MYB axis.")
