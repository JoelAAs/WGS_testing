import numpy as np
import pandas as pd
from py_hla import HLAanalysis

CHROMS = [
    "chr1",
    "chr2",
    "chr3",
    "chr4",
    "chr5",
    "chr6",
    "chr7",
    "chr8",
    "chr9",
    "chr10",
    "chr11",
    "chr12",
    "chr13",
    "chr14",
    "chr15",
    "chr16",
    "chr17",
    "chr18",
    "chr19",
    "chr20",
    "chr21",
    "chr22",
    "chrX",
    "chrY"
]
test = [
    "Hypersensitivity.All~rest",
    "Hypersensitivity.Narrow~rest",
    "Fototoxicitet~rest",
    "SJS.TEN~rest",
    "Cytopenia.All~rest",
    "White.Blood.Cytopenia~rest",
    "Trombocytopeni~rest",
    "Blödning~rest",
    "Cns.Toxicity~rest",
    "Hypo.Siadh~rest",
    "Levertoxicitet~rest",
    "Metabol.rubbning~rest",
    "Njurtoxicitet~rest",
    "Patologisk.fraktur~rest",
    "Pankreatit~rest",
    "Narkolepsi~rest"
]
rule All:
    input:
        #expand("run_folder/no_pc/{assoc}.no_pc.glm.logistic.hybrid", assoc=test),
        expand("run_folder/Plots/gwas/{assoc}_annotated.png", assoc=test)
        #expand("run_folder/Plots/Burden/{assoc}_all.png", assoc=test)
        #"run_folder/all_top.csv",
        #expand("run_folder/Plots/Assocs/{assoc}/gwas.png", assoc= test)
#expand("run_folder/Plots/Chroms/{assoc}/{chrom}_{assoc}_annotated.png", chrom=CHROMS, assoc = test)



rule analyse_hla:
    output:
        csv = "run_folder/HLA/{assoc}.csv"
    run:
        dir = ()
        shell("mkdir -p run_folder/confusion/{wildcards.assoc}")
        hla_filename = "/home/joel/Documents/swedegene/WGS_pipeline/master/HLA-results/analysis_2021-april/input_data/Wadelius_01_v3_2021_01.xlsx"
        pheno_file = "local_data/phenofile"
        covar_file = f"local_data/PCA/{wildcards.assoc}.eigenvec"
        confusion_folder = f"run_folder/confusion/{wildcards.assoc}"
        output_file = f"run_folder/HLA/{wildcards.assoc}.csv"
        exclude_samples = "local_data/exclude_samples.txt"

        hla_obj = HLAanalysis(hla_filename, wildcards.assoc,pheno_file,
            covar_file,confusion_folder,output_file)
        hla_obj.test_all_HLA()
        hla_obj.write_assoc_full()


rule plot_hla:
    input:
        "run_folder/HLA/{assoc}.csv"
    output:
        "run_folder/Plots/HLA/{assoc}.png"
    shell:
        """
        Rscript src/plot_hla.R {input} {wildcards.assoc} {output}
        """


rule plot_structural:
    output:
        "run_folder/Plots/structural/{assoc}.png"
    shell:
        """
        Rscript src/plot_structural.R {wildcards.assoc} {output}
        """


rule collect:
    input:
        snps = "run_folder/Plots/gwas/{assoc}_annotated.png",
        burden = "run_folder/Plots/Burden/{assoc}_all.png",
        hla = "run_folder/Plots/HLA/{assoc}.png",
        struc = "run_folder/Plots/structural/{assoc}.png"
    output:
        snps = "run_folder/Plots/Assocs/{assoc}/gwas.png",
        burden = "run_folder/Plots/Assocs/{assoc}/skato.png",
        hla = "run_folder/Plots/Assocs/{assoc}/HLA.png",
        struc = "run_folder/Plots/Assocs/{assoc}/structural.png"
    shell:
        """
        cp {input.snps} {output.snps}
        cp {input.burden} {output.burden}
        cp {input.hla} {output.hla}
        cp {input.struc} {output.struc}
        """

rule get_all_high_hits:
    input:
        expand("run_folder/high_hits/{assoc}_top.csv", assoc=test)
    output:
        "run_folder/all_top.csv"
    shell:
        """
        head -n 1 {input[0]} > {output[0]}
            for file in {input}; do
                tail -n +2 $file >> {output};
            done 
        """

rule get_tops:
    input:
        "run_folder/exon_distance/concated/{assoc}_annotated.logistic_distance"
    output:
        "run_folder/high_hits/{assoc}_top.csv"
    shell:
        """
        Rscript src/get_high.R {input} {output} {wildcards.assoc}
        """

rule filter_pcs:
    input:
        "local_data/plink/SNP/{assoc}.{assoc}.glm.logistic.hybrid"
    output:
        "run_folder/no_pc/{assoc}.no_pc.glm.logistic.hybrid"
    shell:
        """
        cat {input} | grep -v "PC" > {output}
        """

rule get_variants_per_chrom:
    input:
        variants = "run_folder/no_pc/{assoc}.no_pc.glm.logistic.hybrid"
    output:
        temp(expand("run_folder/chrom_sep/{chrom}_{{assoc}}.logistic", chrom=CHROMS))
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
         annotated = temp("run_folder/exon_distance/{chrom}_{assoc}_annotated.logistic")
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


rule get_phenotable:
    input:
        excluded = "local_data/exclude_samples.txt",
        pheno = "local_data/phenofile"
    output:
        table = "run_folder/phenofile_counts.table"
    run:
        df_pheno = pd.read_csv(input.pheno, sep = "\t")
        df_excluded = pd.read_csv(input.excluded, header=None)
        df_excluded = df_excluded.rename({0:"IID"}, axis = 1)
        df_pheno = df_pheno[~df_pheno["IID"].isin(df_excluded["IID"])]
        with open(output.table, "w") as w:
            pheno_groups = list(df_pheno.columns.values)
            pheno_groups.remove("IID")
            w.write("assoc\tn_cases\tn_controls\n")
            for assoc in pheno_groups:
                w.write(f"{assoc}\t{sum(df_pheno[assoc] == 2)}\t{sum(df_pheno[assoc] == 1)}\n")




#rule filter_variants:
#    input:
#        freq_case = "local_data/freq/{assoc}_2.afreq",
#        freq_control = "local_data/freq/{assoc}_1.afreq",
#        nr_samples = "run_folder/phenofile_counts.table"
#    output:
#        ids = "run_folder/passing_vf/{assoc}_distance_ids"
#    run:
#        with open(input.nr_samples, "r") as f:
#            for line in f:
#                if wildcards.assoc in line:
#                    break
#        _, n_cases, n_controls = line.split("\t")
#        vf_t = [(1/int(n_cases))**(1/2), (1/(975-int(n_controls)))**(1/2)]
#        with open(output.ids,"w") as w:
#            for i in range(2):
#                with open(input[i], "r") as f:
#                    for line in f:
#                        line = line.split("\t")
#                        try:
#                            if float(line[4]) > vf_t[i]:
#                                w.write(f"{line[1]}\n")
#                        except ValueError:
#                            pass



rule filter_variants:
    input:
        freq = "local_data/freq/{assoc}.afreq",
        nr_samples = "run_folder/phenofile_counts.table"
    output:
        ids = "run_folder/passing_vf/{assoc}_distance_ids"
    run:
        with open(input.nr_samples, "r") as f:
            for line in f:
                if wildcards.assoc in line:
                    break
        _, n_cases, n_controls = line.split("\t")
        vf_t = (1/int(n_cases))**(1/2)
        with open(output.ids,"w") as w:
            with open(input.freq, "r") as f:
                for line in f:
                    line = line.split("\t")
                    try:
                        if float(line[4]) > vf_t:
                            w.write(f"{line[1]}\n")
                    except ValueError:
                        pass


rule filter_assocs:
    input:
        passing_variants = "run_folder/passing_vf/{assoc}_distance_ids",
        pvals = "run_folder/exon_distance/concated/{assoc}_annotated_repeats_removed.logistic"
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


rule removeRepeats:
    params:
        repeatmasker = "/home/joel/Documents/Resourses/reapeatMasker_hg19/repeatMasker_sorted_standard.bed"
    input:
        variants = "run_folder/exon_distance/concated/{assoc}_annotated.logistic"
    output:
        variants = temp("run_folder/exon_distance/concated/{assoc}_annotated_repeats_removed.logistic"),
        masked_variants = "run_folder/masked_variants/{assoc}.logistic"
    run:
        order = {str(i):i for i in range(1, 23)}
        order.update({"X":23})
        repeatFile = open(params.repeatmasker, "r")
        logisticFile = open(input.variants, "r")
        passingVariants = open(output.variants, "w")
        maskedVariants = open(output.masked_variants, "w")
        headerRepeat   = repeatFile.readline()
        headerLogistic = logisticFile.readline()

        passingVariants.write(headerLogistic)

        logisticLine = logisticFile.readline()
        repeatLine = repeatFile.readline()
        repChr, repStart, repEnd = repeatLine.split("\t")[0:3]
        logisticChr, logisticPos = logisticLine.split("\t")[0:2]
        repChr = repChr.replace("chr","")

        while repeatLine != "" or logisticLine != "":
            if logisticChr == "Y":
                passingVariants.write(logisticLine)

            elif order[repChr] == order[logisticChr]:
                if int(repEnd) >= int(logisticPos) >= int(repStart):
                    maskedVariants.write(logisticLine)
                    logisticLine = logisticFile.readline()
                    if logisticLine == "":
                        break
                    logisticChr, logisticPos = logisticLine.split("\t")[0:2]

                elif int(repStart) > int(logisticPos):
                    passingVariants.write(logisticLine)
                    logisticLine = logisticFile.readline()
                    if logisticLine == "":
                        break
                    logisticChr, logisticPos = logisticLine.split("\t")[0:2]

                elif int(repStart) < int(logisticPos):
                    repeatLine = repeatFile.readline()
                    if not repeatLine == "":
                        repChr, repStart, repEnd = repeatLine.split("\t")[0:3]
                        repChr = repChr.replace("chr","")

            elif order[repChr] < order[logisticChr]:
                repeatLine = repeatFile.readline()
                if not repeatLine == "":
                    repChr, repStart, repEnd = repeatLine.split("\t")[0:3]
                    repChr = repChr.replace("chr","")

            elif order[repChr] > order[logisticChr]:
                passingVariants.write(logisticLine)
                logisticLine = logisticFile.readline()
                if logisticLine == "":
                    break
                logisticChr, logisticPos = logisticLine.split("\t")[0:2]

        repeatFile.close()
        logisticFile.close()
        passingVariants.close()
        maskedVariants.close()

rule concat_chroms:
    input:
        expand("run_folder/exon_distance/{chrom}_{{assoc}}_annotated.logistic", chrom=CHROMS)
    output:
        temp("run_folder/exon_distance/concated/{assoc}_annotated.logistic")
    shell:
        """
        head -n 1 {input[0]} > {output[0]}
        for file in {input}; do
            tail -n +2 $file >> {output};
        done 
        """

rule plot_burden_single:
    input:
        "local_data/burden_results/burden_dynamic_threshold_full_results.tsv"
    output:
        "run_folder/Plots/Burden/{assoc}_all.png"
    shell:
        """
        Rscript src/plotSingleBurden.R --args {wildcards.assoc}
        """



rule plot_variants:
    input:
        annotated = "run_folder/exon_distance/concated/{assoc}_annotated.logistic_distance"
    output:
        png = "run_folder/Plots/gwas/{assoc}_annotated.png"
    shell:
        """
        Rscript src/nice_manhattan.R --args {input} {output}
        """