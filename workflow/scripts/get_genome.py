#!/usr/bin/python3

# GET GENOME
# -----------------------------------------------------------------------------
#
# This script attempts to download genome sequence (FASTA) and
# genome annotation (GFF / GTF) files from NCBI using the NCBI datasets
# API, or a similar database. Also, the genome sequence is indexed using samtools.
# Alternatively, a FASTA and GFF file can be
# supplied by the user. Input is roughly checked for validity.

from os import path
from BCBio.GFF import GFFExaminer
from BCBio import GFF
from subprocess import getoutput

input_database = snakemake.params["database"]
input_assembly = snakemake.params["assembly"]
input_fasta = snakemake.params["fasta"]
input_gff = snakemake.params["gff"]
output_path = snakemake.output["path"]
output_fasta = snakemake.output["fasta"]
output_gff = snakemake.output["gff"]
output_log = snakemake.log["path"]
log = []
error = []


def check_fasta(input_fasta, log=[], error=[]):
    with open(input_fasta, "r") as fasta_file:
        fasta = fasta_file.read()
    n_items = fasta.count(">")
    if n_items:
        log += [f"Supplied fasta file '{input_fasta}' was found"]
        log += [f"Supplied fasta file contains {n_items} items"]
    else:
        error += ["The supplied fasta file contains no valid entries starting with '>'"]
    return fasta, log, error


def check_gff(input_gff, log=[], error=[]):
    with open(input_gff, "r") as gff_file:
        gff_examiner = GFFExaminer()
        log += [f"Supplied GFF file '{input_gff}' was found"]
        gff_summary = gff_examiner.available_limits(gff_file)
        log += [
            f"Supplied GFF file contains the following items:",
            "--------------------",
        ]
        for item in gff_summary["gff_source_type"]:
            log += ["-".join(item) + " : " + str(gff_summary["gff_source_type"][item])]
    with open(input_gff, "r") as gff_file:
        new_gff = []
        gff_source_type = []
        for i in snakemake.config["get_genome"]["gff_source_type"]:
            gff_source_type += list(i.items())
        limits = dict(gff_source_type=gff_source_type)
        for rec in GFF.parse(gff_file, limit_info=limits):
            for recfeat in rec.features:
                rec_keys = recfeat.qualifiers.keys()
                if not "Name" in rec_keys:
                    if "locus_tag" in rec_keys:
                        recfeat.qualifiers["Name"] = recfeat.qualifiers["locus_tag"]
                    else:
                        error += [
                            "required fields 'Name','locus_tag' missing in *.gff file"
                        ]
                else:
                    if "locus_tag" in rec_keys:
                        recfeat.qualifiers["trivial_name"] = recfeat.qualifiers["Name"]
                        recfeat.qualifiers["Name"] = recfeat.qualifiers["locus_tag"]
                if not "ID" in rec_keys:
                    if "locus_tag" in rec_keys:
                        recfeat.qualifiers["ID"] = recfeat.qualifiers["locus_tag"]
                    elif "Name" in rec_keys:
                        recfeat.qualifiers["ID"] = recfeat.qualifiers["Name"]
                    else:
                        error += [
                            "required fields 'ID','locus_tag' missing in *.gff file"
                        ]
            new_gff += [rec]
    return new_gff, log, error


if input_database.lower() == "ncbi":
    ncbi_result = getoutput(
        f"datasets summary genome accession {input_assembly} --as-json-lines | "
        + "dataformat tsv genome --fields accession,annotinfo-release-date,organism-name"
    )
    if ncbi_result.startswith("Error"):
        error += [ncbi_result]
        error += [
            "The supplied refseq/genbank ID was not valid. Example for correct input: 'GCF_000009045.1'"
        ]
    elif len(ncbi_result) == 0:
        error += [
            "The result from fetching NCBI genome data has zero length. Please check your internet connection!"
        ]
    else:
        ncbi_genome = [
            i.split("\t")
            for i in ncbi_result.split("\n")
            if not (i.startswith("New version") or i.startswith("Warning"))
        ]
        ncbi_genome = dict(zip(ncbi_genome[0], ncbi_genome[1]))
        log += ["Found the following genome(s):"]
        for k in ncbi_genome.keys():
            log += ["{0}: {1}".format(k, ncbi_genome.get(k))]
        refseq_id = ncbi_genome.get("Assembly Accession")
        if not refseq_id.startswith("GCF_"):
            error += ["The RefSeq ID '{0}' has no valid format.".format(refseq_id)]
        ncbi_command = (
            f"datasets download genome accession {refseq_id}"
            + f" --filename {output_path}/database.zip --include genome,gff3; "
            + f"cd {output_path}; unzip database.zip; rm database.zip"
        )
        copy_command = (
            f"cp {output_path}/ncbi_dataset/data/{refseq_id}/*.fna {output_fasta}; "
            + f"cp {output_path}/ncbi_dataset/data/{refseq_id}/genomic.gff {output_gff}"
        )
        index_command = f"samtools faidx {output_fasta}"
        str_out = getoutput(ncbi_command)
        str_cp = getoutput(copy_command)
        # import and check files
        fasta, log, error = check_fasta(output_fasta, log, error)
        index_out = getoutput(index_command)
        gff, log, error = check_gff(output_gff, log, error)
        # write cleaned gff file
        with open(output_gff, "w") as gff_out:
            GFF.write(gff, gff_out)

elif input_database.lower() == "manual":
    if not path.exists(input_fasta):
        error += ["The parameter 'fasta' is not a valid path to a FASTA file"]
    elif not path.exists(input_gff):
        error += ["The parameter 'gff' is not a valid path to a GFF/GTF file"]
    else:
        # import and check files
        fasta, log, error = check_fasta(input_fasta, log, error)
        gff, log, error = check_gff(input_gff, log, error)
        # export fasta and gff files
        with open(output_fasta, "w") as fasta_out:
            fasta_out.write(fasta)
        with open(output_gff, "w") as gff_out:
            GFF.write(gff, gff_out)
else:
    error += ["The parameter 'database' is none of 'ncbi', 'manual'"]

# print error/log messages
if error:
    print("\n".join(error))
    raise ValueError(
        "Location or format of the supplied genome files was not correct, quitting"
    )
else:
    log += ["Module finished successfully\n"]
    log = ["GET_GENOME: " + i for i in log]
    with open(output_log, "w") as log_file:
        log_file.write("\n".join(log))
