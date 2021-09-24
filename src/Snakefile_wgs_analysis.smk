import glob
import pandas as pd
#assocs = [l.replace("data/freq/","").replace("_1.afreq","") for l in glob.glob("data/freq/*_1.afreq")]
assocs = ["CNS.toxicitet~rest_excluded"]
rule all:
    input:
        expand("run_folder/Plots/{assoc}.png", assoc = assocs)

rule filter_pcs:
    input:
        "data/plink/SNP/{assoc}.{assoc}.glm.logistic.hybrid"
    output:
        "run_folder/no_pc/{assoc}.no_pc.glm.logistic.hybrid"
    shell:
        """
        cat {input} | grep -v "PC" > {output}
        """

rule filter_variants:
    input:
        freq_case = "data/freq/{assoc}_2.afreq",
        freq_control = "data/freq/{assoc}_1.afreq",
        nr_samples  = "data/phenofile_counts.table"
    output:
        maf_t = "run_folder/maf_thresholds/{assoc}.csv",
        ids = "run_folder/passing_vf/{assoc}_ids"
    run:
        with open(input.nr_samples, "r") as f:
            for line in f:
                if wildcards.assoc in line:
                    break
        _, n_cases, _, _ = line.split("\t")
        vf_t = [(1/int(n_cases))**(1/2), (1/(975-int(n_cases)))**(1/2)]
        with open(output.maf_t, "w") as w:
            w.write(f"{wildcards.assoc}\t)" + "\t".join([str(l) for l in vf_t]) + "\n")

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
        passing_variants = "run_folder/passing_vf/{assoc}_ids",
        pvals = "run_folder/no_pc/{assoc}.no_pc.glm.logistic.hybrid"
    output:
        variants = "run_folder/excluded/{assoc}.logistic.hybrid"
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
                        w.write(line[1:])


rule get_gene_positions:
    params:
        gene_table = "data/Homo_sapiens.GRCh37.87.genes.tsv"
    input:
        genes = "data/selected_genes/{assoc}.tsv"
    output:
        selected = "run_folder/selected_genes/{assoc}.tsv"
    run:
        gene_definitions = pd.read_csv(params.gene_table, sep="\t")
        selected_genes = pd.read_csv(input.genes, sep="\t")
        gene_pos = gene_definitions[gene_definitions.gene_name.isin(selected_genes.Gene)]
        gene_pos.to_csv(output.selected, index=None, sep="\t")
        if len(gene_pos.index) != len(selected_genes.index):
            raise ValueError("Not all genes were located!")

rule get_variants_in_genes:
    input:
        selected_genes = rules.get_gene_positions.output.selected,
        variants = rules.filter_assocs.output.variants
    output:
        variants_from_selected = "run_folder/selected_genes/variants/{assoc}.logistic.hybrid"
    run:
        gene_def = pd.read_csv(input.selected_genes, sep="\t", dtype={"CHROM":str, "START":int, "END":int})
        chroms =set (gene_def.CHROM)
        with open(output.variants_from_selected,"w") as w:
            with open(input.variants,"r") as f:
                i = 0
                prev_chrom = ""
                current_start = 0
                current_stop = 0

                for linef in f:
                    i += 1
                    if i == 1:
                        w.write(f"{linef}")

                    line = linef.split("\t")
                    if i % 10000 == 0:
                        print(f"{round(i/8845868,ndigits=3)*100} % done")
                    try:
                        chrom = line[0]
                        pos = int(line[1])

                        if chrom != prev_chrom or pos > current_stop:
                            if gene_def.empty:
                                break
                            if chrom in chroms:
                                tmp_gene_def = gene_def[(gene_def.CHROM == chrom)]
                                print(tmp_gene_def.START)
                                tmp_gene_def = tmp_gene_def.sort_values(by=["START"])
                                current_start = list(tmp_gene_def.START)[0]
                                current_stop = list(tmp_gene_def.END)[0]
                                prev_chrom = chrom
                                gene_def = gene_def[gene_def.gene_name != list(tmp_gene_def.gene_name)[0] ]
                                chroms = set(gene_def.CHROM)


                        if current_start <= pos <= current_stop:
                            w.write(f"{linef}")


                    except ValueError:
                        pass



rule plot_stuff:
    input:
        rules.get_variants_in_genes.output.variants_from_selected
    output:
        "run_folder/Plots/{assoc}.png"
    shell:
        """
        Rscript test_manhattan.R --args {input} {output}
        """
