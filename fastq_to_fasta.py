def fastq_to_fasta_efficient(fastq_file, fasta_file, batch_size=1000):
    """Efficient version for processing large files"""
    with open(fastq_file, 'r') as f_in, open(fasta_file, 'w') as f_out:
        batch = []
        seq_count = 0

        for line in f_in:
            batch.append(line.strip())

            # Process one complete sequence record every 4 lines
            if len(batch) == 4:
                # Write FASTA format
                f_out.write(f">{batch[0][1:]}\n")  # Sequence identifier
                f_out.write(f"{batch[1]}\n")  # Sequence

                batch = []  # Clear batch
                seq_count += 1

                # Display progress every certain number of sequences
                if seq_count % batch_size == 0:
                    print(f"Processed {seq_count} sequences...")

        print(f"Conversion completed! Total of {seq_count} sequences processed")


# Usage example
fastq_to_fasta_efficient('', '')