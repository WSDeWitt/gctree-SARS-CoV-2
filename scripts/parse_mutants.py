import pandas as pd
from Bio import AlignIO
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment

ref_seq = AlignIO.read(open(snakemake.input[0]), "fasta")[0]
muts_df = pd.read_csv(snakemake.input[1], index_col='representative_strain',
                      dtype=dict(substitutions=str))
out_aln = MultipleSeqAlignment([])

# we'll use these colors to map location categories
location_colors = {"cat_Wuhan": "red",
                   "cat_other China": "blue",
                   "cat_outside China": "green"}

with (open(snakemake.output[0], "w") as variants_f,
      open(snakemake.output[1], "w") as abundances_f,
      open(snakemake.output[2], "w") as idmap_f,
      open(snakemake.output[3], "w") as locationcts_f):
    for i, idx in enumerate(muts_df.index, 1):
        seq = MutableSeq(ref_seq.seq)
        muts = muts_df.loc[idx, "substitutions"]
        if isinstance(muts, str):
            for mut in muts.split(","):
                ref_base = mut[0]
                pos = int(mut[1:-1]) - 1
                alt_base = mut[-1]
                if seq[pos] != ref_base:
                    raise ValueError(
                            f"variant {mut} is inconsistent "
                            f"with base {ref_base} in reference "
                            "sequence")
                seq[pos] = alt_base
        new_id = "root" if idx == snakemake.params.root.replace("_", "/") else f"seq{i}"
        out_aln.append(SeqRecord(seq, id=new_id, description=""))

        strain_ids = muts_df.loc[idx, "all_strains"].split(", ")
        n_strains = muts_df.loc[idx, "nstrains"]
        assert n_strains == len(strain_ids)
        print(f"{new_id},{n_strains}", file=abundances_f)
        print(f"{new_id},{':'.join(strain_ids)}", file=idmap_f)
        print(f"{new_id}\t" + ",".join(
                  [f"{location_colors[location]}:{muts_df.loc[idx, location]}"
                   for location in location_colors]), file=locationcts_f)

    # remove invariant sites
    variable_sites = [i for i in range(out_aln.get_alignment_length())
                      if len(set(out_aln[:, i])) > 1]
    assert len(variable_sites) > 0
    variable_aln = out_aln[:, variable_sites[0]:variable_sites[0] + 1]
    for site in variable_sites[1:]:
        variable_aln += out_aln[:, site:site + 1]

    # sort so that "root" is first, needed for rerooting
    variable_aln.sort()

    print(format(variable_aln, "phylip"), file=variants_f)
