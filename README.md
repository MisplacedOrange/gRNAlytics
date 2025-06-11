🧬 CRISPR-Cas9 gRNA Specificity Analyzer

This Python tool helps you evaluate and compare the specificity of guide RNAs (gRNAs) designed for CRISPR-Cas9 genome editing. By using the NCBI BLAST database, it identifies potential off-target effects and scores each gRNA based on their match significance, genomic location, and identity to human sequences.

🔍 Features
Accepts multiple user-input gRNA sequences

Runs BLAST searches against the human genome (txid9606)

Scores each gRNA based on:

E-value of off-target hits

Percent identity of matches

Whether matches are on the same chromosome

Transcript/mRNA penalty

Identifies whether genes matched are essential (via AchillesCommonEssentialControls.csv)

Outputs clear results and ranks gRNAs by specificity

ASCII interface and informative printouts

🧑‍🔬 Why Python?
Python provides powerful libraries for:

Bioinformatics – via Biopython

Data handling – with clean parsing and scoring logic

Customization & Accessibility – no expensive lab equipment or microchips needed

This tool can be run on any standard computer using VS Code or any Python IDE.

🧰 Requirements
Install the required packages using pip:

bash
Copy
Edit
pip install biopython numpy
Also make sure you have:

Python 3.7+

Internet access (for BLAST via NCBI)

The file AchillesCommonEssentialControls.csv in the same directory

🚀 How to Run
Clone this repository or download the .py script.

Run the script using:

bash
Copy
Edit
python crispr_grna_analysis.py
Input the number of gRNAs you want to analyze.

Paste each gRNA sequence when prompted.

Review the printed output and results.

📊 Output Example
BLAST results summary for each gRNA

Gene matches and classification (Essential / Non-essential)

Specificity scores out of 100

Final ranking of best-performing gRNA

⚠️ Limitations
This program simulates analysis without a wet lab.

BLAST searches are limited by internet speed and NCBI query limits.

Some gRNAs may return no significant matches, which may appear as errors (but often indicates excellent specificity).


📜 License
This project is for educational and non-commercial research purposes.

