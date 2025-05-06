from Bio import SeqIO

input_fasta = snakemake.input["fasta"]
output_fasta = snakemake.output["fasta"]
log = snakemake.log


def validate_fasta(input_fasta, output_fasta):
    try:
        with open(input_fasta, "r") as fasta_file:
            records = list(SeqIO.parse(fasta_file, "fasta"))
            if not records:
                raise ValueError("FASTA file is empty or improperly formatted.")
            else:
                summary = [f"Validated sequence records for {output_fasta}:"]
                summary += [f"{i.name}: {i.description}" for i in records]
        with open(output_fasta, "w") as validated_file:
            SeqIO.write(records, validated_file, "fasta")
        with open(str(log), "w") as log_file:
            log_file.write("\n".join(summary))
    except Exception as e:
        with open(str(log), "a") as log_file:
            log_file.write(f"Validation failed: {e}\n")
        raise


if __name__ == "__main__":
    validate_fasta(input_fasta, output_fasta)
