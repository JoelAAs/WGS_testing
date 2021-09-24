import numpy as np
import pandas as pd


class GetGenes:
    def __init__(self, gene_tsv):
        self.genes_df = pd.read_csv(gene_tsv, sep="\t")

    @staticmethod
    def _get_first_element_in_list(empty_value, input_list):
        for i, list_value in enumerate(input_list):
            if list_value != empty_value:
                return i, list_value
        return None, None

    def get_distance_of_chrom(self, chrom):
        position_dict = dict()
        chrom_subset_df = self.genes_df[self.genes_df.CHROM == chrom]
        chrom_subset_df = chrom_subset_df.sort_values(by=["START", "END"])
        max_pos = max(chrom_subset_df.END) +1
        chrom_distance = np.full(max_pos, -9, dtype=int)
        chrom_closest_gene = np.empty(max_pos, dtype="O")
        position_dict[chrom] = {
            "pos": chrom_distance,
            "gene": chrom_closest_gene
        }
        for _, (CHROM, START, END, NAME) in chrom_subset_df.iterrows():
            position_dict[chrom]["pos"][START:END+1] = 0
            position_dict[chrom]["gene"][START:END+1] = NAME

        prev_pos = 0
        _, prev_gene = self._get_first_element_in_list("", position_dict[chrom]["gene"])
        first = True
        in_element = False
        for i, distance in enumerate(position_dict[chrom]["pos"]):
            if distance == -9:
                in_element = False
            elif not in_element:
                if first:
                    position_dict[chrom]["gene"][1:i] = prev_gene
                    position_dict[chrom]["pos"][1:i] = range(i - 1, 0, -1)
                    first = False
                else:
                    current_gene = position_dict[chrom]["gene"][i]
                    midpoint = int((i - 1 + prev_pos) / 2)
                    position_dict[chrom]["gene"][prev_pos+1:midpoint+1] = prev_gene
                    position_dict[chrom]["gene"][(midpoint + 1):i] = current_gene

                    position_dict[chrom]["pos"][prev_pos+1:midpoint+1] = range(1, midpoint - prev_pos + 1)
                    position_dict[chrom]["pos"][(midpoint + 1):i] = range(i - midpoint - 1, 0, -1)
                    in_element = True
                    prev_gene = current_gene
            else:
                prev_pos = i

        return position_dict

    def get_some_gene_pickles(self, path):
        for chrom in set(self.genes_df.CHROM):
            chrom_pickle = self.get_distance_of_chrom(chrom)
            np.save(f"{path}/{chrom}.npy", chrom_pickle)