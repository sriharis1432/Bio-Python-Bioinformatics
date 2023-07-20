# Bio-Python-Bioinformatics

# intall packages
!pip install wget
!pip install Bio

# Load the package
import wget

# Replace this URL with the actual file URL you want to download
url = 'https://ftp.ncbi.nlm.nih.gov/genomes/Viruses/Vieuvirus.fn'  
 # Replace this with the desired destination path for the downloaded file.
destination_path = 'Vieuvirus.fn'
# Download the file from the URL to the specified destination path
wget.download(url, destination_path)

print("Download completed!")


# create the bed file with fasta file
from Bio import SeqIO
fasta_file='/content/Vieuvirus.fn'
def fasta_to_bed(fasta_file, bed_file):
    with open(bed_file, 'w') as bed:
        for record in SeqIO.parse(fasta_file, 'fasta'):
            chromosome = record.id
            sequence_length = len(record.seq)
            start_position = 0
            end_position = sequence_length

            bed.write(f"{chromosome}\t{start_position}\t{end_position}\n")

# Replace 'input.fasta' with the path to your FASTA file
# Replace 'output.bed' with the desired output BED file path
fasta_to_bed('/content/Vieuvirus.fn', 'sample.bed')


# load packages
import random
from Bio import SeqIO
from Bio.Seq import Seq

def extract_subsequences(fasta_file, bed_file):
    """
    Extract subsequences from a given primary multi-fasta file based on BED coordinates.

    Parameters:
        fasta_file (str): Path to the primary multi-fasta file.
        bed_file (str): Path to the BED file containing the coordinates for extraction.

    Returns:
        dict: A dictionary containing the extracted subsequences, where keys are scaffold names
              and values are corresponding sequences.
    """
    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
    extracted_subsequences = {}

    with open(bed_file, "r") as bed:
        for line in bed:
            scaffold_name, start_pos, end_pos = line.strip().split()
            start_pos, end_pos = int(start_pos), int(end_pos)
            extracted_subsequences[scaffold_name] = sequences[scaffold_name][start_pos-1:end_pos]

    return extracted_subsequences

def reverse_complement_sequences(extracted_subsequences):
    """
    Reverse complement the extracted subsequences.

    Parameters:
        extracted_subsequences (dict): A dictionary of extracted subsequences.

    Returns:
        dict: A dictionary containing the reverse complemented subsequences,
              where keys are scaffold names and values are corresponding reverse complemented sequences.
    """
    reversed_sequences = {}
    for scaffold_name, sequence in extracted_subsequences.items():
        reversed_sequences[scaffold_name] = sequence.reverse_complement()
    return reversed_sequences

def introduce_random_mutation(reversed_sequences):
    """
    Introduce random 1 nucleotide change in the reversed sequences.

    Parameters:
        reversed_sequences (dict): A dictionary of reverse complemented sequences.

    Returns:
        dict: A dictionary containing the sequences with random one-nucleotide changes,
              where keys are scaffold names and values are corresponding sequences with mutations.
    """
    mutated_sequences = {}
    for scaffold_name, sequence in reversed_sequences.items():
        random_pos = random.randint(0, len(sequence) - 1)
        original_nucleotide = sequence[random_pos]
        mutated_nucleotide = random.choice("ACGT".replace(original_nucleotide, ""))
        mutated_sequence = sequence[:random_pos] + Seq(mutated_nucleotide) + sequence[random_pos + 1:]
        mutated_sequences[scaffold_name] = mutated_sequence
    return mutated_sequences

def insert_processed_subsequences(fasta_file, mutated_sequences, output_file):
    """
    Insert the processed subsequences back into the primary sequence (multifasta).

    Parameters:
        fasta_file (str): Path to the primary multi-fasta file.
        mutated_sequences (dict): A dictionary of sequences with random one-nucleotide changes.
        output_file (str): Path to the output processed multi-fasta file.

    Returns:
        None: The processed multi-fasta file is written to the output_file path.
    """
    with open(output_file, 'w') as output_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            scaffold_name = record.id
            if scaffold_name in mutated_sequences:
                record.seq = mutated_sequences[scaffold_name]
            SeqIO.write(record.seq, output_fasta, "fasta")

# Replace 'Vieuvirus.fn' with the path to your primary multi-fasta file
fasta_file = "Vieuvirus.fn"
# Sample BED file with at least 3 scaffolds and coordinates
bed_file = "sample.bed"
# Output file for the processed multi-fasta
output_file = "processed_" + fasta_file

# Step 1: Extract subsequences based on BED coordinates
extracted_subsequences = extract_subsequences(fasta_file, bed_file)

# Step 2: Reverse complement the extracted subsequences
reversed_sequences = reverse_complement_sequences(extracted_subsequences)

# Step 3: Introduce random nucleotide changes in the reversed sequences
mutated_sequences = introduce_random_mutation(reversed_sequences)

# Step 4: Insert the processed subsequences back into the primary sequence (multifasta)
insert_processed_subsequences(fasta_file, mutated_sequences, output_file)
