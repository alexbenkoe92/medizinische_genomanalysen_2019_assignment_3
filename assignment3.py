#! /usr/bin/env python3

import vcf
import httplib2
from pathlib import Path
import json
from collections import namedtuple

__author__ = 'Alexander Benkö'


##
##
## Aim of this assignment is to annotate the variants with various attributes
## We will use the API provided by "myvariant.info" - more information here: https://docs.myvariant.info
## NOTE NOTE! - check here for hg38 - https://myvariant.info/faq
## 1) Annotate the first 900 variants in the VCF file
## 2) Store the result in a data structure (not in a database)
## 3) Use the data structure to answer the questions
##
## 4) View the VCF in a browser
##

class Assignment3:

    def __init__(self, filename):
        ## Check if pyvcf is installed
        print("PyVCF version: %s" % vcf.VERSION)

        ## Call annotate_vcf_file here
        self.vcf_path = f"{filename}.vcf"

    def annotate_vcf_file(self):
        '''
        - Annotate the VCF file using the following example code (for 1 variant)
        - Iterate of the variants (use first 900)
        - Store the result in a data structure
        :return:
        '''

        ##
        ## Example loop
        ##

        ## Build the connection

        h = httplib2.Http()
        headers = {'content-type': 'application/x-www-form-urlencoded'}

        params_pos = []  # List of variant positions
        with open(self.vcf_path) as my_vcf_fh:
            vcf_reader = vcf.Reader(my_vcf_fh)
            for counter, record in enumerate(vcf_reader):
                params_pos.append(record.CHROM + ":g." + str(record.POS) + record.REF + ">" + str(record.ALT[0]))

                if counter >= 899:
                    break

        ## Build the parameters using the list we just built
        params = 'ids=' + ",".join(params_pos) + '&hg38=true'

        datasource = ["cadd", "dbsnp", "snpeff", "mutdb", "clinvar", "dbsfp"]

        ## Perform annotation
        res, con = h.request('http://myvariant.info/v1/variant', 'POST', params, headers=headers)
        annotation_result = con.decode('utf-8') ## Type String
        self.annotation = json.loads(annotation_result)
        self.annotation_short = []
        for i in self.annotation:
            if not "notfound" in i:
                self.annotation_short.append(i)

        # ## Alternative Way of finding Parameters (without an output file)
        # ## Temporary List, shortened to only include Genes that have been found
        # templist = []
        # for x in self.annotation_short:
        #     for k,v in x.items():
        #         if k in datasource:
        #             templist.append(v)
        #
        # ## GENES & GENENAMES
        # self.gene_names_alt = []
        # for i in templist:
        #     for k,v in i.items():
        #         if type(v) == dict and "genename" in v:
        #             for kx,vx in v.items():
        #                 if kx == "genename":
        #                     self.gene_names_alt.append(vx)
        #         elif type(v) == list:
        #             for vy in v:
        #                 if type(vy) == dict and "genename" in vy:
        #                     self.gene_names_alt.append(vy["genename"])
        # self.gene_names_alt = set(self.gene_names_alt)
        #
        # ## PUTATIVE IMPACT
        # self.imp_modifier_alt = 0
        # for i in templist:
        #     for k,v in i.items():
        #         if type(v) == dict and "putative_impact" in v:
        #             for kx,vx in v.items():
        #                 if kx == "putative_impact" and vx == "MODIFIER":
        #                     self.imp_modifier_alt += 1
        #         elif type(v) == list:
        #             for vy in v:
        #                 if type(vy) == dict and "putative_impact" in vy:
        #                     if vy["putative_impact"] == "MODIFIER":
        #                         self.imp_modifier_alt += 1
        #
        # ## MUTATIONTASTER
        # self.mut_taster_alt = 0
        # for i in templist:
        #     for k,v in i.items():
        #         if type(v) == dict and "mutationtaster" in v:
        #             self.mut_taster_alt += 1
        #         elif type(v) == list:
        #             for vy in v:
        #                 if type(vy) == dict and "mutationtaster" in vy:
        #                     self.mut_taster_alt += 1
        #
        #
        # ## CONSEQUENCE
        # self.cons_nonsyn_alt = 0
        # for i in templist:
        #     for k,v in i.items():
        #         if type(v) == dict and "consequence" in v:
        #             for kx,vx in v.items():
        #                 if kx == "consequence" and vx == "NON_SYNONYMOUS":
        #                     self.cons_nonsyn_alt += 1
        #         elif type(v) == list:
        #             for vy in v:
        #                 if type(vy) == dict and "consequence" in vy:
        #                     if vy["consequence"] == "NON_SYNONYMOUS":
        #                         self.cons_nonsyn_alt += 1
        #
        #
        # print(f"\nGene Names: {self.gene_names_alt}")
        # print(f"Amount of Genes: {len(self.gene_names_alt)}")
        # print(f"\nVariants with Putative Impact = Modifier: {self.imp_modifier_alt}")
        # print(f"\nVariants with Mutationtaster Annotation: {self.mut_taster_alt}")
        # print(f"\nVariants with Consequence = Non-Synonymous: {self.cons_nonsyn_alt}")

        if not Path.cwd().joinpath("annotation.json").exists():
            with open("annotation.json", "w") as af:
                for line in annotation_result:
                    af.write(line)

    def run_analysis(self):
        '''
        Print the name of genes in the annotation data set
        :return:
        '''
        ## Easy Way of finding Parameters
        self.gene_names = []
        self.imp_modifier = 0
        self.mut_taster = 0
        self.cons_nonsyn = 0

        if not Path.cwd().joinpath("annotation.json").exists():
            self.annotate_vcf_file()

        with open("annotation.json", "r") as fh:
            for line in fh:
                if "genename" in line:
                    self.gene_names.append(line)
                if '"putative_impact": "MODIFIER"' in line:
                    self.imp_modifier += 1
                if "mutationtaster" in line:
                    self.mut_taster += 1
                if '"consequence": "NON_SYNONYMOUS"' in line:
                    self.cons_nonsyn += 1

        self.gene_names = set([i.strip().strip(",").strip('"genename": ').strip('"') for i in self.gene_names])
        self.gene_amount = len(self.gene_names)
        print(f"\nGene Names: {self.gene_names}")
        print(f"Amount of Genes: {self.gene_amount}")
        print(f"\nVariants with Putative Impact = Modifier: {self.imp_modifier}")
        print(f"\nVariants with Mutationtaster Annotation: {self.mut_taster}")
        print(f"\nVariants with Consequence = Non-Synonymous: {self.cons_nonsyn}")


    def view_vcf_in_browser(self):
        '''
        - Open a browser and go to https://vcf.iobio.io/
        - Upload the VCF file and investigate the details
        :return:
        '''

        ## Document the final URL here
        print("\nhttps://vcf.iobio.io/?species=Human&build=GRCh38")

    def print_summary(self):
        self.annotate_vcf_file()
        self.run_analysis()


def main():
    print("Assignment 3\n")
    assignment3 = Assignment3("chr16")
    assignment3.print_summary()
    assignment3.view_vcf_in_browser()
    print("\nDone with assignment 3")


if __name__ == '__main__':
    main()
