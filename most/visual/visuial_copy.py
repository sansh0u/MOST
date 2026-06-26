import snapatac2 as snap
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import scanpy as sc

adata = snap.pp.import_fragments(
    fragment_file="/data/Fanyt/ATAC_0209_output/ATAC_sorted.bed.gz",
    chrom_sizes=snap.genome.mm10,
    sorted_by_barcode=False,
)
fig = snap.pl.frag_size_distr(adata, show=False)
fig.update_layout(
    xaxis_title="Fragment size (bp)", yaxis_title="Count",
    width=500,
    height=350,
    margin=dict(
    l=60,
    r=20,
    b=50,
    t=20)
)
fig.update_layout(
    plot_bgcolor='white',  # Remove background color (transparent)
    paper_bgcolor='white',  # Remove the background of the entire figure (transparent)
    xaxis=dict(
        showline=True,  # Show frame on x-axis
        showgrid=False,
        linewidth=0.8,  # Thickness of the frame line
        linecolor='black',  # Color of the frame line
        ticks='outside',  # Place ticks outside the axis line
        tickwidth=1,  # Thickness of the ticks
        tickcolor='black',  # Color of the ticks
        ticklen=4,  # Length of the ticks
        tickmode='auto',  # Can be 'auto' or 'array'
        # Optional: Define specific tick positions
        # tickvals=[0, 1, 2, 3, 4],
        # Optional: Define custom tick labels
        # ticktext=['A', 'B', 'C', 'D', 'E'],
        mirror=True,
        title_font=dict(
            family="Arial",  # Font family
            size=12,  # Font size
            color="black"  # Font color
        ),tickfont=dict(
            family="Arial",  # Font family
            size=10,  # Font size
            color="black"  # Font color
        )
    ),
    yaxis=dict(
        showline=True,  # Show frame on y-axis
        showgrid=False,
        rangemode="tozero",
        tickformat=",",
        linewidth=0.8,  # Thickness of the frame line
        linecolor='black',  # Color of the frame line
        ticks='outside',  # Place ticks outside the axis line
        tickwidth=1,  # Thickness of the ticks
        tickcolor='black',  # Color of the ticks
        ticklen=4,  # Length of the ticks
        tickmode='auto',  # Can be 'auto' or 'array'
        mirror=True,title_font=dict(
            family="Arial",  # Font family
            size=12,  # Font size
            color="black"  # Font color
        ),tickfont=dict(
            family="Arial",  # Font family
            size=10,  # Font size
            color="black"  # Font color
        )
    ),
)
fig.update_xaxes(
    range=[0,700]
)

fig.show()
fig.write_image('fragment_dis.png', scale=6)

snap.metrics.tsse(adata, snap.genome.mm10)

# =========================
# 4. 🔥 TSSe density scatter（论文级替换）
# =========================

x = np.log10(adata.obs["n_fragment"] + 1)
y = adata.obs["tsse"]

# density
xy = np.vstack([x, y])
z = gaussian_kde(xy)(xy)

# sort for clean plot
idx = z.argsort()
x, y, z = x[idx], y[idx], z[idx]

plt.figure(figsize=(5,4))

scat = plt.scatter(
    x, y,
    c=z,
    s=6,
    cmap="turbo",
    edgecolor="none"
)

plt.colorbar(scat, label="Density")

plt.xlabel("log10-transformed unique fragments")
plt.ylabel("TSS enrichment score")

plt.tight_layout()
plt.savefig("tsse_density.png", dpi=300, bbox_inches="tight")
plt.show()


"""
snap.pp.add_tile_matrix(adata)

snap.pp.select_features(adata, n_features=500000, inplace = True)

data2 = adata[:,adata.var['selected']==True]

snap.tl.spectral(data2, n_comps=20)

snap.tl.umap(data2, min_dist=0.02)

sc.settings.set_figure_params(dpi=300, facecolor='white',fontsize=12)
sc.pl.umap(data2, color=["n_fragment",'frac_dup','tsse'], use_raw=False, size=25, wspace=0.25)


snap.pp.knn(data2)
data2

snap.tl.leiden(data2, resolution=1.4, key_added='leiden')
sc.pl.umap(data2,color='leiden', size=25)


## Generate a spatial obsm

data2.uns['spatial']={}
data2.uns['spatial']['whole']={

                           }
data2.obsm['spatial']=pd.DataFrame.to_numpy(data2.obs[['x_coord','y_coord']].astype('int64'))

sc.pl.spatial(data2, color=["leiden"],img_key=None, spot_size=1)

sc.pl.spatial(data2, color=["n_fragment"],img_key=None, spot_size=1, vmax=10000
             )
sc.pl.spatial(data2, color=["tsse"],img_key=None, spot_size=1
             )



palette_cmap=['gainsboro']*10
original_cmap=['steelblue', 'orange', '#279e68', '#d62728', 'darkviolet', '#8c564b', 'navy', '#b5bd61', '#17becf', '#aec7e8']

for i in range(10):
    palette_cmap[i]=original_cmap[i]
    sc.pl.spatial(data2, color=["leiden"],img_key=None, spot_size=0.92, title='Cluster #'+str(i),
             palette=palette_cmap, legend_loc=None)
    palette_cmap=['gainsboro']*10
"""