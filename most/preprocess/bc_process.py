from Bio.SeqIO.QualityIO import FastqGeneralIterator
from gzip import open as gzopen

def atac_bc(input_file,output_file_R1,output_file_R2,seq_start,bc2_start,bc2_end,bc1_start,bc1_end):

    
    with gzopen(input_file, "rt") as in_handle_R1, open(output_file_R1, "w") as out_handle_R1, open(output_file_R2, "w") as out_handle_R2:
        for title, seq, qual in FastqGeneralIterator(in_handle_R1):
            new_seq_R1 = seq[seq_start:]
            new_qual_R1 = qual[seq_start:]
            barcode = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end] # !!! BC2 + BC1
            new_qual_R2 = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end]        
            out_handle_R1.write("@%s\n%s\n+\n%s\n" % (title, new_seq_R1, new_qual_R1))
            out_handle_R2.write("@%s\n%s\n+\n%s\n" % (title, barcode, new_qual_R2))
    

def dbit_bc(input_file,output_file,umi_start,umi_len,bc2_start,bc2_end,bc1_start,bc1_end):


    with gzopen(input_file, "rt") as in_handle:
        with open(output_file, "w") as out_handle:
            for title, seq, qual in FastqGeneralIterator(in_handle):
                new_seq = seq[bc2_start:bc2_end] + seq[bc1_start:bc1_end] + seq[umi_start:umi_start+umi_len]  # BC2 + BC1 + UMI
                new_qual = qual[bc2_start:bc2_end] + qual[bc1_start:bc1_end] + qual[umi_start:umi_start+umi_len]
                out_handle.write("@%s\n%s\n+\n%s\n" % (title, new_seq, new_qual))
    


