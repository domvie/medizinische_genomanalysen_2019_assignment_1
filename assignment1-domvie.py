import mysql.connector
import os
import pybedtools
import pysam


__author__ = 'Dominic Viehboeck'
class Assignment1:
    
    def __init__(self):
        self.gene = "FAM3B"
        self.genome_reference = "hg38"
        #create dictionary out of entry columns
        self.genedict = self.download_gene_coordinates(file_name="test.txt")

        #Set bamfile directories
        self.bamfile = os.path.join(os.getcwd(), "chr21.bam")
        self.baifile = os.path.join(os.getcwd(), "chr21.bam.bai")
        self.samfile = pysam.AlignmentFile(self.bamfile, "rb")
        self.reads = list(self.samfile.fetch("chr21", self.genedict["txStart"], self.genedict["txEnd"]))

    def download_gene_coordinates(self, file_name):
        print("Connecting to UCSC to fetch data")
        
        ## Open connection
        cnx = mysql.connector.connect(host='genome-mysql.cse.ucsc.edu', user='genomep', passwd='password', db=self.genome_reference)
        
        ## Get cursor
        cursor = cnx.cursor()
        
        ## Build query fields
        query_fields = ["refGene.name2",
                        "refGene.name",
                        "refGene.chrom",
                        "refGene.txStart",
                        "refGene.txEnd",
                        "refGene.strand",
                        "refGene.exonCount",
                        "refGene.exonStarts",
                        "refGene.exonEnds"]
        
        ## Build query for gene of interest
        query = 'SELECT DISTINCT %s from refGene' % ','.join(query_fields) + \
        ' WHERE refGene.name2="' + self.gene + '"'

        
        ## Execute query
        cursor.execute(query)

        ## Write to file and create dict with gene entry of interest
        genedict = {}
        with open(file_name, "w") as fh:
            for row in cursor:
                fh.write(str(row) + "\n")
                genedict = {
                    "name2": row[0],
                    "name": row[1],
                    "chrom": row[2],
                    "txStart": row[3],
                    "txEnd": row[4],
                    "strand": row[5],
                    "exonCount": row[6],
                    "exonStarts": row[7],
                    "exonEnds": row[8]
                }

        ## Close cursor & connection
        cursor.close()
        cnx.close()
        fh.close()
        print("Done fetching data")
        return genedict
        
    def get_coordinates_of_gene(self):
        print("Coordinates of gene " + self.gene + ": ")
        print("Start:\t", self.genedict["txStart"], "\nEnd:\t", self.genedict["txEnd"])

    def get_gene_symbol(self):
        print("Gene Symbol: ")
        print(self.genedict["name"])
        print()
                        
    def get_sam_header(self):
        print("Sam Header: \n")
        for key, value in self.samfile.header['HD'].items():
            if key == "SO":
                print("Sorting order of alignments (SO): ", value)
            if key == "VN":
                print("Format version (VN): ", value)
            if key == "GO":
                print("Grouping of alignments (GO): ", value)
                print()

    def get_properly_paired_reads_of_gene(self):
        i = 0
        for read in self.reads:
            if read.is_proper_pair:
                i += 1
        if i == 0:
            print("No properly paired reads found \n")
        else:
            print("Number of properly paired reads: ")
            print(i, "\n")
            
    def get_gene_reads_with_indels(self):
        i = 0
        # Find indels in unmapped reads by using Cigar
        for read in self.reads:
            if not read.is_unmapped:
                cigar = read.cigar
                for (type, length) in cigar:
                    # Insertion or Deletion
                    if (type == 1) or (type == 2):
                        i += 1
        if i == 0:
            print("No gene reads with indels found \n")
        else:
            print("Number of gene reads with indels:")
            print(i,"\n")
        
    def calculate_total_average_coverage(self):
        print("Start calculating total average coverage...\n")
        a = pybedtools.BedTool(self.bamfile)
        b = a.genome_coverage(bg=True)

        i = 0
        average = 0
        for line in b:
            number = float(line[3])
            average += number
            i += 1

        coverage = average/i
        print("Total average coverage:")
        print(coverage, "\n")
        
    def calculate_gene_average_coverage(self):
        print("Start calculating gene average coverage...\n")
        a = pybedtools.BedTool(self.bamfile)
        b = a.genome_coverage(bg=True)

        average = 0
        i = 0

        for line in b:
            number = float(line[3])
            cbeg = int(line[1])
            if cbeg > self.genedict["txStart"]:
                if int(line[2]) <= self.genedict["txEnd"]:
                    average += number
                    i += 1

        coverage = average / i

        print("Total gene average coverage:")
        print(coverage, "\n")
        
    def get_number_mapped_reads(self):
        i = 0
        for read in self.reads:
            if not read.is_unmapped:
                i += 1
        if i == 0:
            print("No mapped reads found \n")
        else:
            print("Number of mapped reads:")
            print(i, "\n")

    def get_region_of_gene(self):
        print("Region of gene:")
        print("Chromosome: ", self.genedict["chrom"])
        print()
        
    def get_number_of_exons(self):
        print("Number of exons: ", self.genedict["exonCount"])
        print()
    
    def print_summary(self):
        self.get_coordinates_of_gene()
        self.get_gene_symbol()
        self.get_sam_header()
        self.get_properly_paired_reads_of_gene()
        self.get_gene_reads_with_indels()
        self.calculate_total_average_coverage()
        self.calculate_gene_average_coverage()
        self.get_number_mapped_reads()
        self.get_region_of_gene()
        self.get_number_of_exons()
        self.samfile.close()

def main():
    print("Assignment 1")
    assignment1 = Assignment1()
    assignment1.print_summary()

    print("Done with assignment 1")
    
        
if __name__ == '__main__':
    print(__author__, "\n")
    main()