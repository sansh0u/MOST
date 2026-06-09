import os
import snapatac2 as snap
def plot_fragment_distribution(cfg):
    output_dir = cfg.out_dir
    adata = snap.pp.import_fragments(
    fragment_file=f"{cfg.out_dir}/{cfg.project}.bed.gz",
    chrom_sizes=getattr(
        snap.genome,
        cfg.reference.genome
    ),
    sorted_by_barcode=False,
)
    fig = snap.pl.frag_size_distr(adata, show=False)
    fig.update_layout(
        xaxis_title="Fragment size (bp)",
        yaxis_title="Count",
        width=500,
        height=350,
        margin=dict(
            l=60,
            r=20,
            b=50,
            t=20
        ),
        plot_bgcolor="white",
        paper_bgcolor="white"
    )

    fig.update_xaxes(
        range=[0,700],
        showline=True,
        showgrid=False,
        linewidth=0.8,
        linecolor="black",
        ticks="outside",
        tickwidth=1,
        ticklen=4,
        title_font=dict(
            family="Arial",
            size=12
        ),
        tickfont=dict(
            family="Arial",
            size=10
        )
    )

    fig.update_yaxes(
        rangemode="tozero",
        tickformat=",",
        showline=True,
        showgrid=False,
        linewidth=0.8,
        linecolor="black",
        ticks="outside",
        tickwidth=1,
        ticklen=4,
        title_font=dict(
            family="Arial",
            size=12
        ),
        tickfont=dict(
            family="Arial",
            size=10
        )
    )
    png_file = os.path.join(
        output_dir,
        "fragment_distribution.png"
    )
    fig.write_image(png_file,scale=6)
