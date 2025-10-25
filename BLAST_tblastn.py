import subprocess
import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
import tempfile
import shutil

# ================================
# Configuration Section
# ================================

# BLAST tool paths
BLAST_PATH = "C:\\Program Files\\NCBI\\blast-2.17.0+\\bin"

# File path configuration
ARABIDOPSIS_PROTEIN_FILE = ""
KANDELIA_GENOME_FILE = ""
OUTPUT_DIR = ""

# BLAST parameters
EVALUE_THRESHOLD = 1e-5
MIN_IDENTITY = 30
TOP_CANDIDATES = 10

# Primer design parameters
PRIMER_LENGTH = 20
PRODUCT_SIZE = 150


# ================================
# Utility Functions
# ================================

def calculate_gc(sequence):
    """Calculate GC content of a sequence"""
    gc_count = sequence.count('G') + sequence.count('C') + sequence.count('g') + sequence.count('c')
    return (gc_count / len(sequence)) * 100 if sequence else 0


def calculate_tm(primer):
    """Calculate primer TM value (simplified method)"""
    if not primer:
        return 0

    a_count = primer.count('A') + primer.count('a')
    t_count = primer.count('T') + primer.count('t')
    g_count = primer.count('G') + primer.count('g')
    c_count = primer.count('C') + primer.count('c')

    return 2 * (a_count + t_count) + 4 * (g_count + c_count)


# ================================
# Core Functions
# ================================

def get_blast_tool(tool_name):
    """Get full path to BLAST tool"""
    if BLAST_PATH and os.path.exists(BLAST_PATH):
        tool_path = os.path.join(BLAST_PATH, f"{tool_name}.exe")
        if os.path.exists(tool_path):
            return tool_path
    return tool_name


def run_blast_analysis(arabidopsis_file, kandelia_file, output_dir):
    """Run complete BLAST analysis pipeline"""
    print("Step 1: Building BLAST database...")

    # Create temporary working directory
    temp_dir = tempfile.mkdtemp()
    db_name = os.path.join(temp_dir, "kandelia_db")

    makeblastdb_cmd = get_blast_tool("makeblastdb")
    cmd = f'"{makeblastdb_cmd}" -in "{kandelia_file}" -dbtype nucl -out "{db_name}"'

    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True, text=True)
        print("Database construction completed")
    except subprocess.CalledProcessError as e:
        print(f"Database construction failed: {e}")
        if hasattr(e, 'stderr') and e.stderr:
            print(f"Error details: {e.stderr}")
        shutil.rmtree(temp_dir)
        return None

    print("Step 2: Running tblastn alignment...")
    blast_output = os.path.join(temp_dir, "tblastn_results.txt")
    tblastn_cmd = get_blast_tool("tblastn")

    cmd = f'"{tblastn_cmd}" -query "{arabidopsis_file}" -db "{db_name}" -out "{blast_output}" -outfmt 6 -evalue {EVALUE_THRESHOLD} -num_threads 4'

    try:
        subprocess.run(cmd, shell=True, check=True)
        print("tblastn alignment completed")
    except subprocess.CalledProcessError as e:
        print(f"tblastn execution failed: {e}")
        shutil.rmtree(temp_dir)
        return None

    # Parse BLAST results
    print("Step 3: Parsing BLAST results...")
    columns = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
               "qstart", "qend", "sstart", "send", "evalue", "bitscore"]

    try:
        # Check if file is empty
        if os.path.getsize(blast_output) == 0:
            print("No BLAST matches found")
            shutil.rmtree(temp_dir)
            return None

        df = pd.read_csv(blast_output, sep='\t', header=None, names=columns)

        # Filter significant matches and sort
        significant_hits = df[
            (df['evalue'] < 1e-5) &
            (df['pident'] > MIN_IDENTITY)
            ].sort_values('evalue').head(TOP_CANDIDATES)

        if significant_hits.empty:
            print("No significant matches found")
            shutil.rmtree(temp_dir)
            return None

        print(f"Found {len(significant_hits)} best candidate genes")

        # Extract candidate gene sequences
        candidate_genes = extract_candidate_genes(significant_hits, kandelia_file)

        # Clean up temporary directory
        shutil.rmtree(temp_dir)

        return significant_hits, candidate_genes

    except Exception as e:
        print(f"Failed to parse BLAST results: {e}")
        shutil.rmtree(temp_dir)
        return None


def extract_candidate_genes(blast_results, genome_file):
    """Extract candidate gene sequences from BLAST results"""
    candidate_genes = []

    try:
        genome_records = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    except Exception as e:
        print(f"Failed to read genome file: {e}")
        return []

    for idx, hit in blast_results.iterrows():
        seq_id = hit['sseqid']
        start = min(hit['sstart'], hit['send'])
        end = max(hit['sstart'], hit['send'])
        strand = '+' if hit['sstart'] < hit['send'] else '-'

        if seq_id in genome_records:
            # Extend region for primer design
            extension = 300
            extended_start = max(1, start - extension)
            extended_end = end + extension

            try:
                sequence = genome_records[seq_id].seq[extended_start - 1:extended_end]

                if strand == '-':
                    sequence = sequence.reverse_complement()
                    # Adjust position information
                    start, end = extended_start, extended_end
                else:
                    start, end = extended_start, extended_end

                candidate_genes.append({
                    'id': f"Candidate_{idx + 1}",
                    'original_location': f"{seq_id}:{hit['sstart']}-{hit['send']}",
                    'extended_location': f"{seq_id}:{start}-{end}",
                    'strand': strand,
                    'evalue': hit['evalue'],
                    'identity': hit['pident'],
                    'sequence': str(sequence),
                    'length': len(sequence),
                    'gc_content': calculate_gc(str(sequence))
                })
            except Exception as e:
                print(f"Failed to extract gene sequence {seq_id}: {e}")
                continue

    return candidate_genes


def design_primers(sequence, gene_id, primer_length=20, product_size=150):
    """Design primers for gene sequence"""
    if len(sequence) < product_size:
        return None

    try:
        # Design forward primer (at sequence start)
        forward_primer = sequence[:primer_length]

        # Design reverse primer (at sequence end, take reverse complement)
        reverse_sequence = sequence[-primer_length:]
        reverse_primer = str(Seq(reverse_sequence).reverse_complement())

        # Calculate product sequence
        product_sequence = sequence[:product_size]

        # Calculate primer properties
        forward_gc = calculate_gc(forward_primer)
        reverse_gc = calculate_gc(reverse_primer)
        product_gc = calculate_gc(product_sequence)

        # Calculate TM values
        forward_tm = calculate_tm(forward_primer)
        reverse_tm = calculate_tm(reverse_primer)

        return {
            'gene_id': gene_id,
            'forward_primer': forward_primer,
            'reverse_primer': reverse_primer,
            'product_size': product_size,
            'forward_tm': forward_tm,
            'reverse_tm': reverse_tm,
            'forward_gc': forward_gc,
            'reverse_gc': reverse_gc,
            'product_gc': product_gc,
            'product_sequence': product_sequence
        }
    except Exception as e:
        print(f"Failed to design primers for gene {gene_id}: {e}")
        return None


def generate_analysis_report(blast_results, candidate_genes, output_dir):
    """Generate analysis report"""
    report_file = os.path.join(output_dir, "analysis_report.txt")

    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("Arabidopsis and Kandelia Homologous Gene Analysis Report\n")
        f.write("=" * 60 + "\n\n")

        f.write("1. Analysis Overview\n")
        f.write("-" * 40 + "\n")
        f.write(f"Arabidopsis query gene: {os.path.basename(ARABIDOPSIS_PROTEIN_FILE)}\n")
        f.write(f"Kandelia genome file: {os.path.basename(KANDELIA_GENOME_FILE)}\n")
        f.write(f"BLAST E-value threshold: {EVALUE_THRESHOLD}\n")
        f.write(f"Minimum identity: {MIN_IDENTITY}%\n")
        f.write(f"Total matches: {len(blast_results)}\n")
        f.write(f"Best candidate genes: {len(candidate_genes)}\n\n")

        f.write("2. Candidate Gene Details\n")
        f.write("-" * 40 + "\n")

        for i, gene in enumerate(candidate_genes, 1):
            f.write(f"\nCandidate Gene #{i}:\n")
            f.write(f"  Gene ID: {gene['id']}\n")
            f.write(f"  Original location: {gene['original_location']}\n")
            f.write(f"  Extended location: {gene['extended_location']}\n")
            f.write(f"  Strand: {gene['strand']}\n")
            f.write(f"  E-value: {gene['evalue']:.2e}\n")
            f.write(f"  Identity: {gene['identity']:.1f}%\n")
            f.write(f"  Sequence length: {gene['length']} bp\n")
            f.write(f"  GC content: {gene['gc_content']:.1f}%\n")

        f.write("\n3. Statistical Analysis\n")
        f.write("-" * 40 + "\n")
        if candidate_genes:
            evalues = [gene['evalue'] for gene in candidate_genes]
            identities = [gene['identity'] for gene in candidate_genes]
            lengths = [gene['length'] for gene in candidate_genes]
            gc_contents = [gene['gc_content'] for gene in candidate_genes]

            f.write(f"Best E-value: {min(evalues):.2e}\n")
            f.write(f"Worst E-value: {max(evalues):.2e}\n")
            f.write(f"Highest identity: {max(identities):.1f}%\n")
            f.write(f"Lowest identity: {min(identities):.1f}%\n")
            f.write(f"Average identity: {sum(identities) / len(identities):.1f}%\n")
            f.write(f"Average sequence length: {sum(lengths) / len(lengths):.0f} bp\n")
            f.write(f"Average GC content: {sum(gc_contents) / len(gc_contents):.1f}%\n")

    print(f"Analysis report generated: {report_file}")
    return report_file


def generate_best_candidates_fasta(candidate_genes, output_dir):
    """Generate FASTA file for best candidate genes"""
    fasta_file = os.path.join(output_dir, "best_candidate_genes.fasta")

    with open(fasta_file, 'w', encoding='utf-8') as f:
        for gene in candidate_genes:
            f.write(f">{gene['id']} | location={gene['extended_location']} | "
                    f"strand={gene['strand']} | evalue={gene['evalue']:.2e} | "
                    f"identity={gene['identity']:.1f}% | length={gene['length']}bp | "
                    f"gc_content={gene['gc_content']:.1f}%\n")
            f.write(f"{gene['sequence']}\n")

    print(f"Best candidate genes FASTA file generated: {fasta_file}")
    return fasta_file


def generate_primer_report(primers, output_dir):
    """Generate primer design report"""
    report_file = os.path.join(output_dir, "primer_design_report.txt")

    with open(report_file, 'w', encoding='utf-8') as f:
        f.write("Kandelia Candidate Genes Primer Design Report\n")
        f.write("=" * 60 + "\n\n")

        f.write("Primer Design Parameters:\n")
        f.write(f"  Primer length: {PRIMER_LENGTH} bp\n")
        f.write(f"  Product size: {PRODUCT_SIZE} bp\n")
        f.write(f"  Designed genes: {len(primers)}\n\n")

        f.write("Primer Design Results:\n")
        f.write("-" * 60 + "\n")

        for i, primer in enumerate(primers, 1):
            f.write(f"\n{i}. Gene: {primer['gene_id']}\n")
            f.write(f"   Forward primer: {primer['forward_primer']}\n")
            f.write(f"     TM value: {primer['forward_tm']:.1f}°C, GC content: {primer['forward_gc']:.1f}%\n")
            f.write(f"   Reverse primer: {primer['reverse_primer']}\n")
            f.write(f"     TM value: {primer['reverse_tm']:.1f}°C, GC content: {primer['reverse_gc']:.1f}%\n")
            f.write(f"   Product size: {primer['product_size']} bp\n")
            f.write(f"   Product GC content: {primer['product_gc']:.1f}%\n")
            f.write(f"   Product sequence: {primer['product_sequence'][:50]}...\n")

        f.write("\nUsage Instructions:\n")
        f.write("-" * 40 + "\n")
        f.write("1. Recommend using professional software to verify primer specificity\n")
        f.write("2. Recommended TM value range: 55-65°C\n")
        f.write("3. Recommended GC content range: 40-60%\n")
        f.write("4. Avoid primer dimers and hairpin structures\n")
        f.write("5. Perform in silico PCR validation before experiments\n")

    print(f"Primer design report generated: {report_file}")
    return report_file


def main():
    """Main function"""
    print("=" * 60)
    print("Arabidopsis and Kandelia Homologous Gene Analysis and Primer Design")
    print("=" * 60)

    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Check input files
    for file_path, desc in [(ARABIDOPSIS_PROTEIN_FILE, "Arabidopsis protein file"),
                            (KANDELIA_GENOME_FILE, "Kandelia genome file")]:
        if not os.path.exists(file_path):
            print(f"Error: Cannot find {desc}: {file_path}")
            return

    # Run BLAST analysis
    blast_results = run_blast_analysis(ARABIDOPSIS_PROTEIN_FILE, KANDELIA_GENOME_FILE, OUTPUT_DIR)
    if not blast_results:
        print("BLAST analysis failed")
        return

    significant_hits, candidate_genes = blast_results

    # Generate analysis report
    analysis_report = generate_analysis_report(significant_hits, candidate_genes, OUTPUT_DIR)

    # Generate best candidate genes FASTA file
    best_candidates = generate_best_candidates_fasta(candidate_genes, OUTPUT_DIR)

    # Design primers
    print("Step 4: Designing primers...")
    primers = []
    for gene in candidate_genes:
        primer = design_primers(gene['sequence'], gene['id'], PRIMER_LENGTH, PRODUCT_SIZE)
        if primer:
            primers.append(primer)

    # Generate primer design report
    primer_report = generate_primer_report(primers, OUTPUT_DIR)

    # Completion summary
    print("\n" + "=" * 60)
    print("Analysis completed!")
    print(f"Output directory: {OUTPUT_DIR}")
    print("Generated files:")
    print(f"  1. {os.path.basename(analysis_report)} - Detailed analysis report")
    print(f"  2. {os.path.basename(best_candidates)} - Candidate gene sequences")
    print(f"  3. {os.path.basename(primer_report)} - Primer design report")
    print("=" * 60)


if __name__ == "__main__":
    # Check required libraries
    try:
        import pandas as pd
        from Bio import SeqIO
        from Bio.Seq import Seq
    except ImportError as e:
        print(f"Error: Missing required Python libraries - {e}")
        print("Please install: pip install biopython pandas")
        exit(1)

    main()