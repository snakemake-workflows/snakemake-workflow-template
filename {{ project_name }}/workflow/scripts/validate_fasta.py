import sys
from Bio import SeqIO

sys.stderr = open(snakemake.log[0], "w", buffering=1)


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
        sys.stderr.write("\n".join(summary))
    except Exception as e:
        sys.stderr.write(f"Validation failed: {e}\n")
        raise


validate_fasta(snakemake.input["fasta"], snakemake.output["fasta"])
