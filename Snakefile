workdir: "build/"

ref_aln_url = "https://raw.githubusercontent.com/jbloom/SARS-CoV-2_PRJNA612766/main/results/ref_genome/ref_genome.fa"
muts_url = "https://raw.githubusercontent.com/jbloom/SARS-CoV-2_PRJNA612766/main/results/phylogenetics/all_alignment.csv"

roots = ["hCoV-19/USA/WA1/2020", "hCoV-19/Guangdong/HKU-SZ-002/2020", "hCoV-19/Shandong/LY005-2/2020"]


rule all:
    input:
        expand("{root}/gctree.log", root=roots)


rule reference:
    """Download reference genome."""
    output:
        "ref_genome.fa"
    shell:
        "curl {ref_aln_url} -O"


rule mutations:
    """Download Jesse's mutation data."""
    output:
        "all_alignment.csv"
    shell:
        "curl {muts_url} -O"


rule variants:
    """Parse mutants csv file to generate an alignment and abundance
    file"""
    input:
        rules.reference.output,
        rules.mutations.output,
    params:
        root = "{root}"
    output:
        "{root}/variants.phy",
        "{root}/abundances.csv",
        "{root}/idmap.csv",
        "{root}/location_color_counts.csv"
    script:
        "scripts/parse_mutants.py"


rule dnapars_config:
    """Make config file to run PHYLIP's dnapars program to infer maximum
    parsimony trees. Use ``--quick`` parameter for less exhaustive tree space
    search.
    """
    input:
        rules.variants.output[0]
    output:
        "{root}/dnapars.cfg"
    shell:
        "mkconfig {input} dnapars --quick > {output}"


rule dnapars:
    """Run dnapars"""
    input:
        rules.dnapars_config.output
    output:
        "{root}/dnapars.log",
        "{root}/outfile",
        "{root}/outtree"
    shell:
        "dnapars < {input} > {output[0]} && mv outfile {output[1]} && mv outtree {output[2]}"


rule gctree:
    """Run gctree to rank maximum parsimony trees using abundance data."""
    input:
        rules.dnapars.output[1],
        rules.variants.output[1],
        rules.variants.output[3]
    output:
        "{root}/gctree.log"
    shell:
        "gctree infer {input[0]} {input[1]} --colormapfile {input[2]} "
        "--outbase {wildcards.root}/gctree.out > {output}"
