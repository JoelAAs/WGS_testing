import numpy as np

CHROMS = [
    "chr10",
    "chr12",
    "chr14",
    "chr16",
    "chr18",
    "chr1",
    "chr21",
    "chr2",
    "chr4",
    "chr6",
    "chr8",
    "chrX",
    "chr11",
    "chr13",
    "chr15",
    "chr17",
    "chr19",
    "chr20",
    "chr22",
    "chr3",
    "chr5",
    "chr7",
    "chr9"
]

test = ["CNS.toxicitet~rest", "Cytopenia.All~rest", "White.Blood.Cytopenia~rest"]
rule All:
    input:
        expand("run_folder/Plots/{assoc}_annotated.png", assoc=test)
        #expand("run_folder/Plots/Chroms/{assoc}/{chrom}_{assoc}_annotated.png", chrom=CHROMS, assoc = test)



rule get_variants_per_chrom:
    input:
        variants = "run_folder/no_pc/{assoc}.no_pc.glm.logistic.hybrid"
    output:
        expand("run_folder/chrom_sep/{chrom}_{{assoc}}.logistic", chrom=CHROMS)
    run:
        chroms_dict = {}
        for chrom in CHROMS:
            chroms_dict.update({
                chrom: open(f"run_folder/chrom_sep/{chrom}_{wildcards.assoc}.logistic", "a")
            })
            chroms_dict[chrom].write("#CHROM\tPOS\tID\tREF\tALT\tA1\tFIRTH?\tTEST\tOBS_CT\tOR\tLOG(OR)_SE\tZ_STAT\tP\n")

        with open(input.variants, "r") as f:
            for line in f:
                c_chr = line.split("\t")[0]
                if c_chr != "#CHROM" and c_chr != "MT" :
                    chroms_dict["chr"+c_chr].write(line)

        for c in chroms_dict:
            chroms_dict[c].close()


rule get_variants_in_genes:
    input:
        distance = "/home/joel/Documents/Resourses/exon_distance_files/{chrom}.npy",
        variants = "run_folder/chrom_sep/{chrom}_{assoc}.logistic"
    output:
         annotated = "run_folder/exon_distance/{chrom}_{assoc}_annotated.logistic"
    run:
        chrom_pos = np.load(input.distance, allow_pickle=True)
        chrom_pos = chrom_pos.item()
        max_len = len(chrom_pos[wildcards.chrom]["pos"])
        with open(input.variants, "r") as f:
            with open(output.annotated, "w") as w:
                first = True
                for line in f:
                    line = line.strip()
                    if first:
                        out_line = line[1:] + "\tDistance\tClosestElement\n"
                        first = False
                    else:
                        line_comp = line.split("\t")
                        pos = int(line_comp[1])
                        if max_len <= pos:
                            dist = pos - max_len
                            gene = chrom_pos[wildcards.chrom]["gene"][-1]
                        else:
                            dist = chrom_pos[wildcards.chrom]["pos"][pos]
                            gene = chrom_pos[wildcards.chrom]["gene"][pos]
                        out_line = line + f"\t{dist}\t{gene}\n"

                    w.write(out_line)



rule filter_variants:
    input:
        freq_case = "data/freq/{assoc}_2.afreq",
        freq_control = "data/freq/{assoc}_1.afreq"
    output:
        ids = "run_folder/passing_vf/{assoc}_distance_ids"
    run:
        vf_t = [0.05, 0.05]
        with open(output.ids,"w") as w:
            for i in range(2):
                with open(input[i], "r") as f:
                    for line in f:
                        line = line.split("\t")
                        try:
                            if float(line[4]) > vf_t[i]:
                                w.write(f"{line[1]}\n")
                        except ValueError:
                            pass

rule filter_assocs:
    input:
        passing_variants = "run_folder/passing_vf/{assoc}_distance_ids",
        pvals = "run_folder/exon_distance/concated/{assoc}_annotated.logistic"
    output:
        variants = "run_folder/exon_distance/concated/{assoc}_annotated.logistic_distance"
    run:
        variants = []
        with open(input.passing_variants, "r") as f:
            variants += [l.strip() for l in f.readlines()]
        variants = set(variants)

        with open(output.variants, "w") as w:
            with open(input.pvals, "r") as f:
                i = 0
                for line in f:
                    i += 1
                    line_s = line.split("\t")
                    if line_s[2] in variants:
                        variants.remove(line_s[2])
                        w.write(line)
                    if i == 1:
                        w.write(line)

rule concat_chroms:
    input:
        expand("run_folder/exon_distance/{chrom}_{assoc}_annotated.logistic", chrom=CHROMS, assoc = test)
    output:
        "run_folder/exon_distance/concated/{assoc}_annotated.logistic"
    shell:
        """
        head -n 1 {input[0]} > {output[0]}
        for file in {input}; do
            tail -n +2 $file >> {output};
        done 
        """

rule plot_variants:
    input:
        annotated = "run_folder/exon_distance/concated/{assoc}_annotated.logistic_distance"
    output:
        png = "run_folder/Plots/{assoc}_annotated.png"
    shell:
        """
        Rscript src/nice_manhattan.R --args {input} {output}
        """