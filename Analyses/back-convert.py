import sys

# Made w/ help from Gemini

# --- CONFIGURATION ---
mapping_file = "coordinate_map.csv"
input_vcf = "LYJA.superContigs.calledGTs.SNPs.filtered.missing.cp.vcf"
output_vcf = "LYJA.superContigs.calledGTs.SNPs.filtered.missing.cp.back-convert.vcf"

# 1. Load the map into a searchable dictionary of lists
# Structure: { 'pseudo_1': [(start, end, contig_name), ...], ... }
coords = {}
print("Loading mapping file...")
with open(mapping_file, 'r') as f:
    next(f) # skip header
    for line in f:
        pseudo, p_start, p_end, orig_id, orig_len = line.strip().split(',')
        if pseudo not in coords:
            coords[pseudo] = []
        coords[pseudo].append((int(p_start), int(p_end), orig_id))

# 2. Process the VCF
print("Converting VCF coordinates...")
with open(input_vcf, 'r') as vcf_in, open(output_vcf, 'w') as vcf_out:
    for line in vcf_in:
        # Pass headers through unchanged
        if line.startswith("#"):
            vcf_out.write(line)
            continue
        
        cols = line.split('\t')
        chrom = cols[0]
        pos = int(cols[1])
        
        # Find the original contig
        found = False
        if chrom in coords:
            # Look for the range that contains 'pos'
            # (Note: For 8m contigs, you could use binary search here for more speed)
            for start, end, orig_id in coords[chrom]:
                if start <= pos <= end:
                    new_pos = pos - start + 1
                    cols[0] = orig_id
                    cols[1] = str(new_pos)
                    vcf_out.write('\t'.join(cols))
                    found = True
                    break
        
        if not found:
            # This happens if a SNP falls in the 'N' gap buffer
            # Usually these are filtered out or ignored
            continue

print(f"Done! Saved to {output_vcf}")
