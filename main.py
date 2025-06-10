from Bio.Blast import NCBIWWW, NCBIXML
import math
import numpy as np
import difflib
import re
from typing import List, Tuple

# Display ASCII art
print("""


██████╗░██╗░██████╗░█████╗░██╗░░░░░░█████╗░██╗███╗░░░███╗███████╗██████╗░██╗
██╔══██╗██║██╔════╝██╔══██╗██║░░░░░██╔══██╗██║████╗░████║██╔════╝██╔══██╗╚═╝
██║░░██║██║╚█████╗░██║░░╚═╝██║░░░░░███████║██║██╔████╔██║█████╗░░██████╔╝░░░
██║░░██║██║░╚═══██╗██║░░██╗██║░░░░░██╔══██║██║██║╚██╔╝██║██╔══╝░░██╔══██╗░░░
██████╔╝██║██████╔╝╚█████╔╝███████╗██║░░██║██║██║░╚═╝░██║███████╗██║░░██║██╗
╚═════╝░╚═╝╚═════╝░░╚════╝░╚══════╝╚═╝░░╚═╝╚═╝░░░░░╚═╝╚══════╝╚═╝░░╚═╝╚═╝



█▄░█ █▀█ ▀█▀   ▄█ █▀█ █▀█ ▀░▄▀   ▄▀█ █▀▀ █▀▀ █░█ █▀█ ▄▀█ ▀█▀ █▀▀
█░▀█ █▄█ ░█░   ░█ █▄█ █▄█ ▄▀░▄   █▀█ █▄▄ █▄▄ █▄█ █▀▄ █▀█ ░█░ ██▄


""")

# Get number of gRNAs from user
numberof_grnas = int(input("How many gRNAs are you comparing? \n").strip())

grna_sequences = []

# Get all gRNA sequences from user
for i in range(numberof_grnas):
    seq = input(f"Enter sequence for gRNA #{i + 1}: \n").strip().upper()
    grna_sequences.append(seq)

print("Starting BLAST search...")

# Load essential genes from CSV file
essential_genes = set()
with open("AchillesCommonEssentialControls.csv", "rt") as file:
    for line in file:
        gene_raw = line.strip().split(",")[0]
        gene_clean = gene_raw.split("(")[0].strip().upper()
        essential_genes.add(gene_clean)

def extract_gene_info(title: str) -> Tuple[str, str]:
    # Extract gene name and symbol from BLAST title with multiple patterns
    gene_patterns = [
        r"gene[:\s]+([A-Z0-9\-_]+)",      #google was used to help write this portion cuz it was rlly hard and we couldn't figure it out 
        r"\(([A-Z0-9\-_]+)\)",
        r"symbol[:\s]+([A-Z0-9\-_]+)",
        r"LOC\d+[,\s]+([A-Z0-9\-_]+)"
    ]
    #google help stop
    for pattern in gene_patterns:
        match = re.search(pattern, title, re.IGNORECASE)
        if match:
            gene_candidate = match.group(1).upper()
            # Check if gene is in essential genes set #refernced the ahillescommonsessentials set
            if gene_candidate in essential_genes:
                return gene_candidate, "Essential"
            else:
                # Try fuzzy matching for partial matches or any direct matches too i think idrmb tbh
                close_matches = difflib.get_close_matches(
                    gene_candidate, essential_genes, n=1, cutoff=0.8
                )
                if close_matches:
                    return close_matches[0], "Essential (fuzzy match)"
                else:
                    return gene_candidate, "Non-essential/Not_Known"
    
    return "Not_Known", "Not_Known"
#if gene is unknown say its unknown
def calculate_specificity_score(alignments: List) -> int:
    # Calculte the scores 
    base_score = 100
    penalty = 0
    
    if not alignments:
        return base_score  # No off-targets found is good
    
    for i, alignment in enumerate(alignments[:10]):  # Analyze top 10 hits
        hsp = alignment.hsps[0]
        
        # Penalty based on e-value (lower e-value means higher chance of cutting off target so higher penalty
        if hsp.expect < 1e-10:
            penalty += 20
        elif hsp.expect < 1e-5:
            penalty += 10
        elif hsp.expect < 0.01:
            penalty += 5
        else:
            penalty += 1
        
        # Penalty based on sequence identity percentage
        identity_percent = (hsp.identities / hsp.align_length) * 100
        if identity_percent > 90:
            penalty += 15
        elif identity_percent > 80:
            penalty += 10
        elif identity_percent > 70:
            penalty += 5
        
        # Additional penalty for transcript matches (potential off-targets)
        if "mRNA" in alignment.title or "transcript" in alignment.title:
            penalty += 3
        
        # Reduce penalty weight for later matches (less important)
        penalty = penalty * (0.8 ** i)
    
    return max(0, base_score - int(penalty))

scores = []

# Process each gRNA sequence
for i, sequence in enumerate(grna_sequences):
    print(f"\nRunning BLAST for gRNA #{i + 1}...\n")
    
    # Run BLAST with enhanced parameters for better analysis
    result_handle = NCBIWWW.qblast(
        program="blastn-short",
        database="nt",
        sequence=sequence,
        entrez_query="txid9606[ORGN]",  # Human sequences only
        hitlist_size=50,  # Get more results for comprehensive analysis
        word_size=7  # Smaller word size better for short gRNA sequences
    )

    # Save BLAST results to XML file
    with open(f"blast_result_{i + 1}.xml", "w") as out_file:
        out_file.write(result_handle.read())

    # Parse BLAST results
    with open(f"blast_result_{i + 1}.xml") as result_file:
        blast_record = NCBIXML.read(result_file)

    # Calculate enhanced specificity score
    specificity_score = calculate_specificity_score(blast_record.alignments)

    if not blast_record.alignments:
        print("No matches found - excellent specificity!")
    else:
        print("Top matches:")
        # Analyze top matches with enhanced gene detection
        for j, alignment in enumerate(blast_record.alignments[:10]):
            hsp = alignment.hsps[0]
            title = alignment.title
            
            # Determine match type
            if "mRNA" in title or "transcript" in title:
                match_type = "Transcript/mRNA"
            else:
                match_type = "Genomic"

            # Enhanced e-value assessment
            if hsp.expect < 1e-10:
                evalue_status = "❗ Very significant"
                evalue_display = f"{hsp.expect:.2e}"
            elif hsp.expect < 1e-5:
                evalue_status = "⚠️  Significant"
                evalue_display = f"{hsp.expect:.2e}"
            elif hsp.expect < 0.01:
                evalue_status = "❓ Moderate"
                evalue_display = f"{hsp.expect:.2e}"
            else:
                evalue_status = "✅ Low significance"
                evalue_display = f"{hsp.expect:.2e}"

            # Extract gene information with enhanced detection
            gene_name, gene_status = extract_gene_info(title)
            
            # Calculate identity percentage
            identity_percent = (hsp.identities / hsp.align_length) * 100

            # Display match information
            print(f"\n--- Match #{j + 1} ---")
            print(f"Type: {match_type}")
            print(f"Gene: {gene_name} ({gene_status})")
            print(f"Title: {title[:80]}...")
            print(f"Length: {alignment.length} bp")
            print(f"E-value: {evalue_display} ({evalue_status})")
            print(f"Identity: {hsp.identities}/{hsp.align_length} ({identity_percent:.1f}%)")
            print(f"Match sequence: {hsp.sbjct[:60]}...")

    # Store the calculated specificity score
    scores.append(specificity_score)
    print(f"\nSpecificity score for gRNA #{i + 1}: {specificity_score}/100")

# Determine the best gRNA based on specificity scores
best_index = scores.index(max(scores))
print("\n" + "="*60)
print(f"🏆 BEST gRNA: gRNA #{best_index + 1}")
print(f"   Sequence: {grna_sequences[best_index]}")
print(f"   Specificity Score: {scores[best_index]}/100")
print("="*60 + "\n")

# Display all scores for comparison
print("All gRNA specificity scores:")
for i, score in enumerate(scores):
    print(f"gRNA #{i + 1}: {score}/100")