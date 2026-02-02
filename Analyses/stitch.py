from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

# Made w/ help from Gemini

# --- CONFIGURATION ---
input_fasta = "LYJA_v0.1.scafSeq.fasta"
output_fasta = "LYJA.stitched_reference.fasta"
mapping_file = "coordinate_map.csv"

min_contig_len = 100        # Filter out very small fragments
target_pseudo_size = 100_000_000  # 100 MB target size
gap_size = 100               # Number of Ns between contigs
# ----------------------

def create_pseudochromosomes():
    with open(mapping_file, "w") as map_out:
        # Header for the lookup table
        map_out.write("pseudo_chr,pseudo_start,pseudo_end,original_contig,original_len\n")
        
        pseudo_count = 1
        current_pseudo_seq = ""
        current_pseudo_len = 0
        
        # Generator to handle 8 million records without crashing RAM
        records = SeqIO.parse(input_fasta, "fasta")
        
        stitched_records = []

        for record in records:
            orig_len = len(record.seq)
            if orig_len < min_contig_len:
                continue
            
            # If current pseudochromosome is full, save it and start a new one
            if current_pseudo_len + orig_len + gap_size > target_pseudo_size:
                new_rec = SeqRecord(Seq(current_pseudo_seq), id=f"pseudo_{pseudo_count}", description="")
                yield new_rec
                
                pseudo_count += 1
                current_pseudo_seq = ""
                current_pseudo_len = 0
            
            # Calculate coordinates
            start_pos = current_pseudo_len + 1
            end_pos = current_pseudo_len + orig_len
            
            # Update mapping file
            map_out.write(f"pseudo_{pseudo_count},{start_pos},{end_pos},{record.id},{orig_len}\n")
            
            # Append sequence and gap
            current_pseudo_seq += str(record.seq) + ("N" * gap_size)
            current_pseudo_len += orig_len + gap_size

        # Don't forget the last batch!
        if current_pseudo_seq:
            yield SeqRecord(Seq(current_pseudo_seq), id=f"pseudo_{pseudo_count}", description="")

# Write the new FASTA file
SeqIO.write(create_pseudochromosomes(), output_fasta, "fasta")
print("Done! Pseudochromosomes and mapping file generated.")
