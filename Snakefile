workdir: "build/"

root = "hCoV-19/USA/WA1/2020"
# root = "hCoV-19/Guangdong/HKU-SZ-002/2020"
# root = "hCoV-19/Shandong/LY005-2/2020"

rule all:
    input:
        "variants.fa",


rule reference:
    output:
        "ref_genome.fa"
    shell:
        "curl https://raw.githubusercontent.com/jbloom/SARS-CoV-2_PRJNA612766"
        "/main/results/ref_genome/ref_genome.fa -O"

rule mutations:
    output:
        "all_alignment_no_filter_rare.csv"
    shell:
        "curl https://raw.githubusercontent.com/jbloom/SARS-CoV-2_PRJNA612766"
        "/main/results/phylogenetics/all_alignment_no_filter_rare.csv -O"

rule variants:
    input:
        rules.reference.output,
        rules.mutations.output,
    output:
        "variants.fa",
    run:
        import pandas as pd
        from Bio import AlignIO
        from Bio.SeqRecord import SeqRecord
        from Bio.Align import MultipleSeqAlignment

        ref_seq = AlignIO.read(open(input[0]), "fasta")[0]
        muts_df = pd.read_csv(input[1], index_col='representative_strain',
                              dtype=dict(substitutions=str))
        out_aln = MultipleSeqAlignment([])

        with open(output[0], "w") as variants_f:

            for i, idx in enumerate(muts_df.index):
                seq = ref_seq.seq
                muts = muts_df.loc[idx, "substitutions"]
                strain_ids = muts_df.loc[idx, "all_strains"].split(", ")
                assert muts_df.loc[idx, "nstrains"] == len(strain_ids)
                if isinstance(muts, str):
                    for mut in muts.split(","):
                        ref_base = mut[0]
                        pos = int(mut[1:-1]) - 1
                        alt_base = mut[-1]
                        if seq[pos] == ref_base:
                            seq[pos] == alt_base
                        else:
                            raise ValueError(
                                    f"variant {mut} is inconsistent "
                                    f"with base {ref_base} in reference "
                                    "sequence")

                new_id = f"seq{i}"
                out_aln.append(SeqRecord(seq, id=new_id, description=idx))
                if idx == root:
                    out_aln.append(SeqRecord(seq, id="root", description=""))

            print(format(out_aln, "fasta"), file=variants_f)
