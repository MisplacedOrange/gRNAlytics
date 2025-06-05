from Bio.Blast import NCBIWWW, NCBIXML
import random
import numpy
import math



# Ask how many gRNAs you want to compare
print("""
╔═╗░╔╦═══╦════╗░░░╔╗╔═══╦═══╦═══╦═══╦═══╦═══╦═══╦═╗░╔╦════╗  ╔═══╦═══╦═══╦╗░╔╦═══╦═══╦════╦═══╗
║║╚╗║║╔═╗║╔╗╔╗║░░╔╝║║╔═╗║╔═╗║╔═╗║╔══╣╔═╗║╔═╗║╔══╣║╚╗║║╔╗╔╗║  ║╔═╗║╔═╗║╔═╗║║░║║╔═╗║╔═╗║╔╗╔╗║╔══╝
║╔╗╚╝║║░║╠╝║║╚╝░░╚╗║║║║║║║║║║╚═╝║╚══╣╚═╝║║░╚╣╚══╣╔╗╚╝╠╝║║╚╝  ║║░║║║░╚╣║░╚╣║░║║╚═╝║║░║╠╝║║╚╣╚══╗
║║╚╗║║║░║║░║║░░░░░║║║║║║║║║║║╔══╣╔══╣╔╗╔╣║░╔╣╔══╣║╚╗║║░║║░░  ║╚═╝║║░╔╣║░╔╣║░║║╔╗╔╣╚═╝║░║║░║╔══╝
║║░║║║╚═╝║░║║░░ ╔╝╚╣╚═╝║╚═╝║║░░║╚══╣║║╚╣╚═╝║╚══╣║░║║║░║║░░  ║╔═╗║╚═╝║╚═╝║╚═╝║║║╚╣╔═╗║░║║░║╚══╗
╚╝░╚═╩═══╝░╚╝░░  ╚══╩═══╩═══╩╝░░╚═══╩╝╚═╩═══╩═══╩╝░╚═╝░╚╝░░  ╚╝░╚╩═══╩═══╩═══╩╝╚═╩╝░╚╝░╚╝░╚═══╝
""")
numberof_grnas = int(input("How many gRNAs are you comparing? \n").strip())

grna_sequences = []

# Collect the sequences based on the number specified
for i in range(numberof_grnas):
    seq = input(f"Enter sequence for gRNA #{i + 1}: \n").strip().upper()
    grna_sequences.append(seq)

print("Starting BLAST search...")

# Run BLAST for each gRNA sequence in the list
for i, sequence in enumerate(grna_sequences):
    print(f"\nRunning BLAST for gRNA #{i + 1}...\n")
    result_handle = NCBIWWW.qblast(
        program="blastn",
        database="nt",
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

            # Simple check to guess if it is transcript or genomic DNA
            if "mRNA" in title or "transcript" in title:
                match_type = "Transcript/mRNA"
            else:
                match_type = "Genomic"
            # Simple check to flag e-value
            if hsp.expect < 0.01:
                evalue = "❗ " + f"{hsp.expect:.2e}"
            else:
                evalue = "❓ " + f"{hsp.expect:.2e}"

            #formats and parses the info  
            print("\n--- Match ---")
            print("Type:", match_type)
            print("Title:", title)
            print("Length:", alignment.length)
            print("E-value:", evalue)               
            print("Match Snippet:", hsp.sbjct[:60] + "...")

