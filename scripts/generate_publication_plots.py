import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

base_dir = "/media/ramesh/LGI/myb-01072026/NLR–NAC–MYB-Axis"
plot_dir = os.path.join(base_dir, "publication_plots")
table_dir = os.path.join(base_dir, "publication_tables")
script_dir = os.path.join(base_dir, "publication_scripts")

os.makedirs(plot_dir, exist_ok=True)
os.makedirs(table_dir, exist_ok=True)
os.makedirs(script_dir, exist_ok=True)

res_dir = os.path.join(base_dir, "analysis_results")

# Load Datasets
df_meta = pd.read_csv(os.path.join(res_dir, "master_target_genes_meta.csv"))
df_modules = pd.read_csv(os.path.join(res_dir, "wgcna_coexpression_modules_310_genes.csv"))
df_motifs = pd.read_csv(os.path.join(res_dir, "promoter_cis_regulatory_motifs_2kb.csv"))
df_expr = pd.read_csv(os.path.join(res_dir, "unified_310_genes_expression_matrix.csv"), index_col=0)

# Set global publication plot style
plt.style.use('seaborn-v0_8-whitegrid')
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10
plt.rcParams['axes.titlesize'] = 12
plt.rcParams['axes.labelsize'] = 11

# ==========================================
# Figure 1: Gene Family Composition & Abundance
# ==========================================
fig, ax = plt.subplots(figsize=(6, 4.5), dpi=300)
fam_counts = df_meta['Family'].value_counts()
colors = ['#1f77b4', '#2ca02c', '#ff7f0e'] # MYB, NLR, NAC colors

bars = ax.bar(fam_counts.index, fam_counts.values, color=colors, edgecolor='black', linewidth=1.2, width=0.55)

for bar in bars:
    yval = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2.0, yval + 2, f"{int(yval)}", ha='center', va='bottom', fontweight='bold', fontsize=11)

ax.set_ylabel('Number of Identified Genes', fontweight='bold')
ax.set_title('Figure 1. Gene Family Distribution (NLR–NAC–MYB Axis)', fontweight='bold', pad=12)
ax.set_ylim(0, max(fam_counts.values) + 20)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig1_png = os.path.join(plot_dir, "Figure_1_Gene_Family_Distribution.png")
fig1_pdf = os.path.join(plot_dir, "Figure_1_Gene_Family_Distribution.pdf")
plt.savefig(fig1_png, dpi=300)
plt.savefig(fig1_pdf)
plt.close()
print("Saved Figure 1.")

# ==========================================
# Figure 2: WGCNA Co-expression Module Composition
# ==========================================
fig, ax = plt.subplots(figsize=(8, 5), dpi=300)
mod_ct = pd.crosstab(df_modules['Module_Name'], df_modules['Family'])

mod_ct.plot(kind='bar', stacked=True, ax=ax, color=['#1f77b4', '#ff7f0e', '#2ca02c'], edgecolor='black', linewidth=0.8)

ax.set_xlabel('WGCNA Co-expression Module', fontweight='bold', labelpad=10)
ax.set_ylabel('Number of Co-expressed Genes', fontweight='bold')
ax.set_title('Figure 2. WGCNA Co-expression Module Composition Across Families', fontweight='bold', pad=15)
ax.legend(title='Gene Family', frameon=True, facecolor='white')
plt.xticks(rotation=15, ha='right')

for p in ax.patches:
    width, height = p.get_width(), p.get_height()
    if height > 5:
        x, y = p.get_xy() 
        ax.text(x + width/2, y + height/2, f"{int(height)}", ha='center', va='center', color='white', fontweight='bold', fontsize=9)

plt.tight_layout()
fig2_png = os.path.join(plot_dir, "Figure_2_WGCNA_Module_Composition.png")
fig2_pdf = os.path.join(plot_dir, "Figure_2_WGCNA_Module_Composition.pdf")
plt.savefig(fig2_png, dpi=300)
plt.savefig(fig2_pdf)
plt.close()
print("Saved Figure 2.")

# ==========================================
# Figure 3: Cis-Regulatory Promoter Syntax Heatmap
# ==========================================
fig, ax = plt.subplots(figsize=(7, 4.5), dpi=300)
motif_means = df_motifs.groupby('Family')[['NAC_binding_site', 'MYB_binding_site', 'W_box_Immune', 'ABRE_Drought']].mean()
motif_means.columns = ['NAC Binding Site', 'MYB Binding Site', 'W-Box (Immune)', 'ABRE (ABA/Drought)']

sns.heatmap(motif_means, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={'label': 'Mean Motif Density / 2.0 kb Promoter'}, ax=ax, linewidths=1.0, linecolor='gray')

ax.set_title('Figure 3. Mean Cis-Regulatory Motif Density (2.0 kb Promoters)', fontweight='bold', pad=15)
ax.set_ylabel('Gene Family', fontweight='bold')

plt.tight_layout()
fig3_png = os.path.join(plot_dir, "Figure_3_Promoter_Motif_Heatmap.png")
fig3_pdf = os.path.join(plot_dir, "Figure_3_Promoter_Motif_Heatmap.pdf")
plt.savefig(fig3_png, dpi=300)
plt.savefig(fig3_pdf)
plt.close()
print("Saved Figure 3.")

# ==========================================
# Figure 4: Multi-Family Stress Expression Heatmap
# ==========================================
fig, ax = plt.subplots(figsize=(8, 6), dpi=300)
# Top 10 Hub genes per family
top_hubs = df_modules.groupby('Family').apply(lambda x: x.nlargest(10, 'kME_Hub_Score'))['Gene_ID'].tolist()

expr_sub = df_expr.loc[df_expr.index.isin(top_hubs)]

sns.heatmap(expr_sub, cmap='vlag', center=0, cbar_kws={'label': 'Expression Z-Score'}, ax=ax, yticklabels=True)
ax.set_title('Figure 4. Expression Profiles of Top NLR, NAC, and MYB Hub Genes', fontweight='bold', pad=15)
ax.set_ylabel('Key Hub Genes', fontweight='bold')
ax.tick_params(axis='y', labelsize=7)

plt.tight_layout()
fig4_png = os.path.join(plot_dir, "Figure_4_Hub_Genes_Expression_Heatmap.png")
fig4_pdf = os.path.join(plot_dir, "Figure_4_Hub_Genes_Expression_Heatmap.pdf")
plt.savefig(fig4_png, dpi=300)
plt.savefig(fig4_pdf)
plt.close()
print("Saved Figure 4.")

# Copy tables to publication_tables
df_meta.to_csv(os.path.join(table_dir, "Table_S1_Target_Genes_Metadata.csv"), index=False)
df_modules.to_csv(os.path.join(table_dir, "Table_S2_WGCNA_Coexpression_Modules.csv"), index=False)
df_motifs.to_csv(os.path.join(table_dir, "Table_S3_Promoter_Motif_Density_2kb.csv"), index=False)
print("Copied all publication tables.")
