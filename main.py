from Bio.Blast import NCBIWWW, NCBIXML
import random
import numpy
import math
import difflib  # Added for fuzzy matching will hopefully be more accurate with the gene flagging

# Ask how many gRNAs I want to compare
print("""




██████╗░██╗░██████╗░█████╗░██╗░░░░░░█████╗░██╗███╗░░░███╗███████╗██████╗░██╗
██╔══██╗██║██╔════╝██╔══██╗██║░░░░░██╔══██╗██║████╗░████║██╔════╝██╔══██╗╚═╝
██║░░██║██║╚█████╗░██║░░╚═╝██║░░░░░███████║██║██╔████╔██║█████╗░░██████╔╝░░░
██║░░██║██║░╚═══██╗██║░░██╗██║░░░░░██╔══██║██║██║╚██╔╝██║██╔══╝░░██╔══██╗░░░
██████╔╝██║██████╔╝╚█████╔╝███████╗██║░░██║██║██║░╚═╝░██║███████╗██║░░██║██╗
╚═════╝░╚═╝╚═════╝░░╚════╝░╚══════╝╚═╝░░╚═╝╚═╝╚═╝░░░░░╚═╝╚══════╝╚═╝░░╚═╝╚



█▄░█ █▀█ ▀█▀   ▄█ █▀█ █▀█ ▀░▄▀   ▄▀█ █▀▀ █▀▀ █░█ █▀█ ▄▀█ ▀█▀ █▀▀
█░▀█ █▄█ ░█░   ░█ █▄█ █▄█ ▄▀░▄   █▀█ █▄▄ █▄▄ █▄█ █▀▄ █▀█ ░█░ ██▄


""")
numberof_grnas = int(input("How many gRNAs are you comparing? \n").strip())

grna_sequences = []

# Collect the sequences based on the number specified
for i in range(numberof_grnas):
    seq = input(f"Enter sequence for gRNA #{i + 1}: \n").strip().upper()
    grna_sequences.append(seq)

print("Starting BLAST search...")

essential_genes = set()
with open("AchillesCommonEssentialControls.csv", "rt") as file:
    for line in file:
        gene_raw = line.strip().split(",")[0]
        gene_clean = gene_raw.split("(")[0].strip().upper()  # removes stuff like " (6059)"
        essential_genes.add(gene_clean)

# Run BLAST for each gRNA sequence in the list
for i, sequence in enumerate(grna_sequences):
    print(f"\nRunning BLAST for gRNA #{i + 1}...\n")
    result_handle = NCBIWWW.qblast(
        program="blastn",
        database="refseq_genomic",
        sequence=sequence,
        entrez_query="txid9606[ORGN]"
    )

    # Save the result with a unique filename to avoid overwriting
    with open(f"blast_result_{i + 1}.xml", "w") as out_file:
        out_file.write(result_handle.read())

    # Parse the BLAST results from the saved XML file
    with open(f"blast_result_{i + 1}.xml") as result_file:
        blast_record = NCBIXML.read(result_file)

    if not blast_record.alignments:
        print("No matches found.")
    else:
        print("Top matches:")
        for alignment in blast_record.alignments[:10]:
            hsp = alignment.hsps[0]
            title = alignment.title

            if "mRNA" in title or "transcript" in title:
                match_type = "Transcript/mRNA"
            else:
                match_type = "Genomic"
            
            # flagging poor matches

            if hsp.expect < 0.01:
                evalue = "❗ " + f"{hsp.expect:.2e}"
            else:
                evalue = "❓ " + f"{hsp.expect:.2e}"

            # Fuzzy gene matching
            title_upper = title.upper()
            best_match = difflib.get_close_matches(
                title_upper, essential_genes, n=1, cutoff=0.6
            )

            if best_match:
                gene_status = "✅ Essential gene (Good to target)"
            else:
                gene_status = "⚠️ Non-essential or unknown"

            # Prints out Data Table

            print("\n--- Match ---")
            print("Type:", match_type)
            print("Title:", title)
            print("Length:", alignment.length)
            print("E-value:", evalue)               
            print("Match Snippet:", hsp.sbjct[:60] + "...")
            print("Gene status:", gene_status)
