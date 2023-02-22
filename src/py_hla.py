import warnings

warnings.filterwarnings('ignore')

import pandas                   as pd
import numpy                    as np
import statsmodels.api          as sm
import statsmodels.formula.api  as smf
from statsmodels.tools.sm_exceptions import PerfectSeparationError
import scipy.stats              as stats
from collections import Counter
from pandas.core.common import flatten
from sys import stderr


class HLAanalysis:
    def __init__(
            self, hla_filename, assoc, pheno_file, covar_file, confusion_folder, output_file, digits=4,
            covar_col=""
    ):

        self.digits = digits
        self.confusion_folder = confusion_folder
        self.output_file = output_file

        self.phenotype_df = pd.read_csv(pheno_file, sep="\t")
        self.associations = [assoc, ]

        self.covar_df = pd.read_csv(covar_file, sep="\t")
        #self.covar_df = self.covar_df.drop(['FID'], axis=1)
        self.covar_df = self.covar_df.set_index('IID')
        if covar_col == "":
            self.covar_columns = ["PC1", "PC2", "PC3", "PC4"]
        elif covar_col is False:
            self.covar_columns = []
        else:
            self.covar_columns = covar_col

        self.hla_df = self.read_and_reformat_hla(hla_filename)
        self.merged_data = self.merge_phenodata()
        self.genetic_model = None
        self.adr_subset = None
        self.adr = None
        self.hlagroup = None
        self.allele = None

        self.full_output = None

    # ######### READ AND FORMAT #########
    def shorten(self, geno):
        # removes genotypes shorter than 4 digits, and those ending in G
        dig = int(self.digits / 2)

        if len(geno) == 2:
            new = []
            for gen in geno:
                if gen == gen:
                    cols = gen.strip().split(':')
                    if len(cols) < 2 or cols[-1].endswith('G'):
                        return [np.nan, np.nan]
                    else:
                        new.append(':'.join(cols[0:dig]))
                else:
                    return [np.nan, np.nan]
            return new
        else:
            return [np.nan, np.nan]

    def read_and_reformat_hla(self, hlafile):
        # read in excel file with hla results and reformat to extend sample lines if several alleles for one sample

        hla = pd.read_excel(hlafile, dtype='object', index_col=0, engine="openpyxl")

        old = list(hla.index)
        new = []
        i = 0
        for entry in old:
            if type(entry) is float:
                new.append(new[i - 1])
            else:
                new.append(entry)
            i += 1

        hla.insert(0, 'SampleName', new)

        col = [c for c in hla.columns if c.startswith('HLA')]
        new_hla = hla.groupby('SampleName')[col].agg(list).reset_index()
        new_hla.set_index('SampleName', inplace=True)

        hla = new_hla.applymap(self.shorten)
        return hla

    def exclude_samples(self, sample_file):
        with open(sample_file, "r") as f:
            exclude_samples = [l.strip for l in f.readline()]
        self.hla_df = self.hla_df.drop(exclude_samples)


    def get_hla_dict(self):
        # get a dict of all unique HLA types present in the data

        hla_dict = self.hla_df.to_dict()
        for entry in hla_dict:
            s = set(list(flatten(hla_dict[entry].values())))
            t = set(filter(lambda x: x == x and type(x) != int, s))
            hla_dict[entry] = t

        return hla_dict

    def merge_phenodata(self):
        # read in phenotype data and PC data and merge with hla data

        merged_phe = pd.merge(self.hla_df, self.phenotype_df, how='inner', left_on='SampleName', right_on='IID')
        merged_phe.set_index('IID')

        merged = pd.merge(merged_phe, self.covar_df, how='inner', left_on='IID', right_on='IID')
        merged.set_index('IID')

        return merged

    # ####### Build genetic model #####
    def set_adr_subset(self, adr):
        # removes all samples designated -9 for case control status, meaning they should not be in the analysis
        self.adr_subset = self.merged_data[self.merged_data[adr] != -9]

    @staticmethod
    def find_additive(geno, allele):
        if pd.isnull(geno).any():  # sets the model value to nan if there are missing hla alleles
            return np.nan
        return geno.count(allele)

    @staticmethod
    def find_dominant(geno, allele):
        if pd.isnull(geno).any():
            return np.nan
        g = 1 if allele in geno else 0
        return g

    @staticmethod
    def find_recessive(geno, allele):
        if pd.isnull(geno).any():
            return np.nan
        c = geno.count(allele)
        return 1 if c == 2 else 0

    def build_model(self, adr, hlagroup, allele, genetic_model):
        self.adr = adr
        self.hlagroup = hlagroup
        self.allele = allele
        if genetic_model not in ['additive', 'allelic', "dominant", "recessive"]:
            raise ImportError(f"Genetic model model {genetic_model} not supported!")
        self.genetic_model = genetic_model

        # assign each genotype a 0,1,or 2 depending on which model to test
        self.set_adr_subset(adr)  # first remove any samples that has a -9 in the ADR column

        if self.genetic_model in ['additive', 'allelic']:
            self.adr_subset['model'] = self.adr_subset[self.hlagroup].apply(self.find_additive, allele=allele)
            self.adr_subset.dropna(inplace=True)

        elif self.genetic_model == 'dominant':
            self.adr_subset['model'] = self.adr_subset[self.hlagroup].apply(self.find_dominant, allele=allele)
            self.adr_subset.dropna(inplace=True)

        elif self.genetic_model == 'recessive':
            self.adr_subset['model'] = self.adr_subset[self.hlagroup].apply(self.find_recessive, allele=allele)
            self.adr_subset.dropna(inplace=True)

    def rebuild_model(self, adr=None, hlagroup=None, allele=None, genetic_model=None):
        self.adr = (adr if adr is not None else self.adr)
        self.hlagroup = (hlagroup if hlagroup is not None else self.hlagroup)
        self.allele = (allele if allele is not None else self.allele)
        self.genetic_model = (genetic_model if genetic_model is not None else self.genetic_model)

        self.build_model(self.adr, self.hlagroup, self.allele, self.genetic_model)

    # ##### Single test #####
    def assoc_test(self):
        # restructure input data and run stat tests

        x = self.adr_subset['model'].values
        y = self.adr_subset[self.adr].values
        z = Counter(zip(x, y))

        header = "{}:{}:{}".format(self.hlagroup, self.allele, self.adr)
        confusion_filename = f"{self.confusion_folder}/{header}.{self.genetic_model}.txt"
        if len(set(x)) < 2 or len(set(y)) < 2:
            with open(confusion_filename, "w+") as out:
                out.write('No association results, not enough data\n\n')
            oddr, pval = "NA"

        elif self.genetic_model == 'additive':
            self.write_output_additive(confusion_filename, z)
            oddr, pval = self.stat_cochran_armitage()

        elif self.genetic_model == 'allelic':
            w = self.count_alleles(z)
            self.write_output_rest(confusion_filename, w)
            oddr, pval = self.stat_fisher(w)
        else:
            self.write_output_rest(confusion_filename, z)
            oddr, pval = self.stat_logistic()

        return {
            "HLA_GROUP": self.hlagroup, "ALLELE": self.allele,
            "TEST": self.adr, "GENETIC_MODEL": self.genetic_model,
            "ODDR": oddr, "PVALUE": pval
        }

    def stat_logistic(self):
        # http://www.statsmodels.org/devel/glm.html, http://www.science.smith.edu/~jcrouser/SDS293/labs/lab4-py.html

        di = {1: 0, 2: 1}
        self.adr_subset['converted'] = self.adr_subset[self.adr].map(di)  # remap case control status to 0 or 1
        # self.adr_subset['model'] = sm.add_constant(self.adr_subset['model'])

        formula = ["converted ~ model"]
        formula += self.covar_columns
        formula = "+".join(formula)
        glm_model = smf.glm(formula=formula, data=self.adr_subset, family=sm.families.Binomial())
        try:
            result = glm_model.fit()  # model='nm'
        except PerfectSeparationError:
            stderr.write(
                f"Perfect seperation error for logistic regression on {self.allele}: {self.adr}."
                f" Genetic model {self.genetic_model}\n")
            return "NA", "NA"
        # TODO: GET ODDS
        return "-", result.pvalues["model"]

    def stat_cochran_armitage(self):
        # https://www.statsmodels.org/dev/contingency_tables.html

        tab = pd.crosstab(self.adr_subset[self.adr], self.adr_subset['model'])
        tab = tab.loc[:, [2, 1, 0]]
        table = sm.stats.Table(tab)
        rslt = table.test_nominal_association()
        # cum_or = table.cumulative_oddsratios

        return "-", rslt.pvalue

    @staticmethod
    def count_alleles(z):
        # count alleles from the table with genotype counts. So two alleles per sample

        w = Counter()
        w.update({(1, 1): z[(2, 1)] * 2 + z[(1, 1)]})
        w.update({(0, 1): z[(1, 1)] + z[(0, 1)] * 2})
        w.update({(1, 2): z[(2, 2)] * 2 + z[(1, 2)]})
        w.update({(0, 2): z[(1, 2)] + z[(0, 2)] * 2})
        return w

    @staticmethod
    def stat_fisher(z):
        # calculate p value
        oddsratio, pvalue = stats.fisher_exact([[z[(1, 2)], z[(0, 2)]], [z[(1, 1)], z[(0, 1)]]])
        return oddsratio, pvalue

    # ##### Multiple tests #####
    def test_all_HLA(self, models=None):
        if models is None:
            models = ['additive', 'allelic', "dominant", "recessive"]

        result_list = []
        hla_dict = self.get_hla_dict()
        for hlagroup in hla_dict:
            for allele in hla_dict[hlagroup]:
                stderr.write(f"Running tests for {allele} variants\n")
                for adr in self.associations:
                    for genetic_model in models:
                        self.build_model(
                            adr=adr, hlagroup=hlagroup, allele=allele, genetic_model=genetic_model
                        )
                        result_list.append(self.assoc_test())
        self.full_output = pd.DataFrame(result_list)

    # ##### Format output #####
    @staticmethod
    def fix(a, c, z):
        # format table more nicely with percentage

        m = max([u[0] for u in z]) + 1
        total = 0
        i = 0
        while i < m:
            total += z[(i, c)]
            i += 1
        total = 0.1 if total == 0 else total
        per = round(a / total * 100)
        per = '0%' if per > 100 else str(per) + '%'

        s = '{:} ({:})'.format(a, per)
        return s

    def write_output_additive(self, outfilename, z):
        # build the tables for the output
        with open(outfilename, "w+") as out:
            out.write("{}\t{}\t{}\n".format(self.hlagroup, self.allele, self.adr))
            out.write('|{:-^10}|{:-^20}|{:-^20}|{:-^20}|\n'.format(
                '', self.allele + ' x2', self.allele + ' x1', 'Other HLA alleles'))
            out.write('|{:10}|{:^20}|{:^20}|{:^20}|\n'.format(
                'Controls', self.fix(z[(2, 1)], 1, z), self.fix(z[(1, 1)], 1, z), self.fix(z[(0, 1)], 1, z)))
            out.write('|{:10}|{:^20}|{:^20}|{:^20}|\n'.format(
                'Cases', self.fix(z[(2, 2)], 2, z), self.fix(z[(1, 2)], 2, z), self.fix(z[(0, 2)], 2, z)))

    def write_output_rest(self, outfilename, z):
        with open(outfilename, "w+") as out:
            out.write("{}\t{}\t{}\n".format(self.hlagroup, self.allele, self.adr))
            out.write('|{:-^10}|{:-^20}|{:-^20}|\n'.format('', self.allele, 'Other HLA alleles'))
            out.write(
                '|{:10}|{:^20}|{:^20}|\n'.format('Controls', self.fix(z[(1, 1)], 1, z), self.fix(z[(0, 1)], 1, z)))
            out.write('|{:10}|{:^20}|{:^20}|\n'.format('Cases', self.fix(z[(1, 2)], 2, z), self.fix(z[(0, 2)], 2, z)))

    def write_assoc_full(self):
        if self.full_output is None:
            raise ValueError("No output generated, run 'test_all_HLA()'")
        else:
            self.full_output.to_csv(self.output_file, sep="\t", index=False)


