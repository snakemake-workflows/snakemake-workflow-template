rule copy_file:
    input:
        config["input"]
    output:
        "results/output.txt"
    shell:
        "cp {input} {output}"
