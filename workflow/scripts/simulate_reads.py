from os import path
from snakemake.shell import shell

# define default options
read_length = snakemake.params.get("read_length", 70)
read_number = snakemake.params.get("read_number", 1e5)
mutations = snakemake.params.get("mutations", 0.001)
random_reads = snakemake.params.get("random_reads", 0.05)
output_type = snakemake.params.get("output_type", 1)
output_prefix = path.commonprefix(snakemake.output).rstrip(".bwa.read")
extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# get input files
fasta = snakemake.input.get("fasta")

# run tool
shell(
    "dwgsim "
    " -1 {read_length}"  # length simulated read 1
    " -2 {read_length}"  # length simulated read 2
    " -N {read_number}"  # total number reads per sample
    " -r {mutations}"  # frequency of mutations
    " -y {random_reads}"  # frequency of random reads
    " -o {output_type}"  # output type, '1' for separate files per read
    " {extra}"  # additional options
    " {fasta}"  # input fasta file
    " {output_prefix}"  # output prefix
    " {log}"
)
