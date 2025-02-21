import os
import numpy as np
import matplotlib.pyplot as plt

def fastq_read_num(fp_fastq):
    """Calculate entropy for a sequence of numbers
 
    Args:
        sequence : contig sequence
 
    Returns:
        float: entropy
 
    """
    sample_names = []
    read_num_fq=[]


    for fp in fp_fastq:

        folder_path = os.path.dirname(fp)
        folder = os.path.basename(folder_path)
        sample_names.append(folder)
        
        with open(fp, 'r') as f:

            count=0
            for line in f:
                if line.startswith("@"):
                    count+=1
                
            read_num_fq.append(count)



    return read_num_fq,sample_names

def render_mpl_table(data, col_width=3.0, row_height=0.625, font_size=12,
                     header_color='#ABC9EA', row_colors=['w'], edge_color='w',
                     bbox=[0, 0, 1, 1], header_columns=0,
                     ax=None, **kwargs):
    if ax is None:
        size = (np.array(data.shape[::-1]) + np.array([0, 1])) * np.array([col_width, row_height])
        fig, ax = plt.subplots(figsize=size)
        ax.axis('off')
    mpl_table = ax.table(cellText=data.values, bbox=bbox, colLabels=data.columns,cellLoc="center", **kwargs)
    mpl_table.auto_set_font_size(False)
    mpl_table.set_fontsize(font_size)

    for k, cell in mpl_table._cells.items():
        cell.set_edgecolor(edge_color)
        if k[0] == 0 or k[1] < header_columns:
            cell.set_text_props(weight='bold', color='w')
            cell.set_facecolor(header_color)
        else:
            cell.set_facecolor(row_colors[k[0]%len(row_colors) ])
    return ax.get_figure(), ax
