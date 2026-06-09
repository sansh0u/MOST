
import snapatac2 as snap

adata = snap.pp.import_fragments(
    fragment_file="/mnt/d/prog/HYS_sample1_sorted.bed.gz",
    chrom_sizes=snap.genome.hg38,
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
