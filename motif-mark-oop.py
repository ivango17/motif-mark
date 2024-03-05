#!/usr/bin/env python

# Author: Ian VanGordon
# Date: 02/18/2024

'''
This program is used to visualize introns/exons and protien binding motifs in mRNA transcripts or DNA.
'''

import argparse
import cairo
import math
import re

def get_args():
    parser = argparse.ArgumentParser(description="This program is used to visualize protien binding motifs in sequences.")
    parser.add_argument("-f", "--fasta", help="Where is the file containing the FASTA sequence?", type=str)
    parser.add_argument("-m", "--motifs", help="What is the destination for the output?", type=str)
    return parser.parse_args()


##################################################################################################################################################################################################################################
# Nucleotide dictionary
##################################################################################################################################################################################################################################
'''
This dictionary uses ambig nucleotides as keys and returns the regular expression of all the bases that the key can be.
'''

ambig_nucleotides = {
    "N" : "[NAGCTU]",
    "R" : "[RAG]",
    "Y" : "[YTCU]",
    "K" : "[KGTU]",
    "M" : "[MAC]",
    "S" : "[SGC]",
    "W" : "[WATU]",
    "B" : "[BCGTU]",
    "D" : "[DAGTU]",
    "H" : "[HACG]",
    "V" : "[VACG]",
    "U" : "[TU]",
    "T" : "[TU]"
}


##################################################################################################################################################################################################################################
# Functions
##################################################################################################################################################################################################################################
'''
Functions to parse input files and create object lists.
'''


def parse_fasta(filepath):
    '''Parses a FASTA file and creates a list of transcript objects from header and sequence.'''
    transcripts_list = []

    with open(filepath) as fa_file:
        readname = ""
        seq = ""
        for line in fa_file:
            line = line.strip()
            if line.startswith(">") and readname == "":
                line = line.split()
                readname = line[0][1:]
            elif line.startswith(">"):
                transcripts_list.append(Transcript(readname, seq))
                seq = ""
                line = line.split()
                readname = line[0][1:]
            else:
                seq += line
        transcripts_list.append(Transcript(readname, seq))

    return transcripts_list


def parse_motif_txt(filepath):
    '''Parses a text file of motifs with only one motif on each line of the file. Returns a list of motif objects.'''
    motifs_list = []

    with open(args.motifs) as motif_file:
        for line in motif_file:
            line = line.strip()
            if line != "":
                motifs_list.append(Motif(line.upper()))

    return motifs_list


##################################################################################################################################################################################################################################
# Classes
##################################################################################################################################################################################################################################
'''
Here the following classes are defined: Motif, Transcript, MotifsImage.
'''

class Motif:
    '''This object will hold information on motif seqs, regex for motif searches, and motif lengths.'''

    def __init__(self, motif_seq):
        self.seq = motif_seq
        self.length = len(motif_seq)
        self.re_seq = self.modular_motif()

    def __repr__(self):
        '''Returns motif sequence as representation of motif object'''
        return f"motif_{self.seq}"

    def modular_motif(self):
        '''Takes a list of motifs and returns a dictionary of ReGexs for motifs as keys that contain ambiguous nucleotides and unmodified motifs that do not contain ambiguous nucleotides. The values are the original motif.'''
        new_motif = ""
        for i in range(len(self.seq)):
            if self.seq[i] in ambig_nucleotides.keys():
                new_motif += ambig_nucleotides[self.seq[i]]
            else:
                new_motif += self.seq[i]
        return new_motif


class Transcript:
    '''This object will hold information about an mRNA read including: sequence, read name, exon locations, and motif locations. '''

    def __init__(self, readname, sequence):
        self.readname = readname
        self.sequence = sequence
        self.exons = [(match.start() + 1) for match in re.finditer(r'[A-Z]', self.sequence)]                                                                                           # Adjusted for "true" position in sequence
        self.length = len(self.sequence)

    def __repr__(self):
        '''Returns read name as representation of transcript object'''
        return f"transcript_{self.readname}"

    def find_motifs(self, motifs_list):
        '''Returns a dictionary (motif : indices) of indices where motifs are present.'''
        self.motif_locations = {}
        for motif in motifs_list:
            self.motif_locations[motif] = [list(range((match.start() + 1), (match.start() + motif.length + 1))) for match in re.finditer('(?={0})'.format(motif.re_seq), self.sequence.upper())] # Adjusted for "true" position in sequence
            self.motif_locations[motif] = [x for l in self.motif_locations[motif] for x in l]                                                                                           # Unlisting the list of lists for each motif

        self.position_motif = {}
        for position in range(1, (len(self.sequence) + 1)):
            self.position_motif[position] = []
            for motif in self.motif_locations.keys():
                if (position in self.motif_locations[motif]) and (motif not in self.position_motif[position]):
                    self.position_motif[position].append(motif.seq)


class MotifsImage:
    '''This object will hold attributes and methods needed to generate a visual representaton introns/exons and protien motif binding sites in sequences.'''

    def __init__(self, name, transcripts_list, motifs_list):

        # Setting image name and assigning transcripts and motifs objects list
        self.name = name
        self.transcripts = transcripts_list
        self.motifs = motifs_list

        # Creating parameters for image dimensions and spacing
        self.scale = 2
        self.num_images = len(transcripts_list)
        self.title_len = 25 + (20 * len(self.motifs))
        self.longest_len = max(transcripts_list, key=lambda x: x.length).length
        self.foot_len = 120
        self.spacing = 200
        self.margins = 50

        # Image dimensions
        self.width = (self.margins * 2) + (self.longest_len * self.scale) 
        self.height = (self.num_images * self.spacing) + self.title_len + self.foot_len

        # Axis attributes
        self.tick_spacing = (self.longest_len * self.scale) // 4
        self.ticklabels = ["0", str(round(self.longest_len / 4)), str(round((self.longest_len / 2))), str(round((self.longest_len * 3) / 4)), str(self.longest_len)]
        self.tickoffset = [0, 12, 12, 12, 15]

        # Generating color pallete
        self.colors = [
            (0, 0.5, 0),        # Green
            (1, 0, 0),          # Red
            (0, 0, 1),          # Blue
            (0.5, 0, 0.5),      # Purple
            (0, 0.75, 1),       # Teal
            (1, 0, 1),          # Pink
            (0.6, 0.2, 0.4),    # Maroon
            (1, 0.4, 0),        # Orange
            (1, 1, 0),          # Yellow
            (1, 0, 0.5)         # Salmon
            ]

        # Assigning colors to motifs
        self.color_pairs = {}
        for co_count, motif in enumerate(motifs_list):
            self.color_pairs[motif.seq] = self.colors[co_count]

    def __repr__(self):
        '''Returns the name of the image as a representation of the image object.'''
        return f"image_{self.name}"

    def create_image(self):
        '''Initializes the image with dimensions from attributes.'''
        self.surface = cairo.SVGSurface(self.name, self.width, self.height)
        self.context = cairo.Context(self.surface)

        self.context.rectangle(0, 0, self.width, self.height)
        self.context.set_source_rgba (255, 255, 255)
        self.context.fill()

    def add_seqlines(self):
        '''Adds sequence lines depending on length of transcript.'''
        self.context.set_line_width(5)
        self.context.set_source_rgba (0, 0, 0)
        for i in range(len(self.transcripts)):
            self.context.move_to(self.margins, (self.title_len + 150 + (self.spacing * i)))
            self.context.line_to((self.scale * self.transcripts[i].length + self.margins), (self.title_len + 150 + (self.spacing * i)))
            self.context.stroke()

    def add_exons(self):
        '''Adds black boxes where exons are present.'''
        self.context.set_line_width(2)
        self.context.set_source_rgb(0, 0, 0)
        for cur_tran, tran in enumerate(self.transcripts):
            for i in range(len(tran.exons)):
                self.context.move_to((self.margins + (self.scale * tran.exons[i])), (self.title_len + 130 + (self.spacing * cur_tran)))
                self.context.line_to((self.margins + (self.scale * tran.exons[i])), (self.title_len + 150 + (self.spacing * cur_tran)))
                self.context.stroke()

    def add_motifs(self):
        '''Adds color coordinated bars representing where motifs are present.'''
        for cur_tran, tran in enumerate(self.transcripts):
            for pos in tran.position_motif.keys():
                if tran.position_motif[pos]:
                    mot_prop = 30 / len(tran.position_motif[pos])
                    for mot_count, motif in enumerate(tran.position_motif[pos]):
                        self.context.set_line_width(2)
                        self.context.move_to((self.margins + (self.scale * pos)), (self.title_len + (150 - mot_prop) + (mot_prop * (mot_count + 1)) + (self.spacing * cur_tran)))
                        self.context.line_to((self.margins + (self.scale * pos)), (self.title_len + 150 + (mot_prop * (mot_count + 1)) + (self.spacing * cur_tran)))
                        self.context.set_source_rgba(*self.color_pairs[tran.position_motif[pos][mot_count]])
                        self.context.stroke()

    def add_readlabels(self):
        '''Adds sequence name to the image.'''
        for tran_count, tran in enumerate(self.transcripts):
            self.context.set_source_rgb(0, 0, 0)
            self.context.set_font_size(30)
            self.context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            self.context.move_to(self.margins, (self.title_len + 110 + (self.spacing * tran_count)))
            self.context.show_text(f"{tran.readname} ({tran.length} bases)")

    def add_ticks(self):
        '''Adds x-axis ticks to show length of reads.'''
        for i in range(len(self.ticklabels)):
            self.context.set_line_width(3)
            self.context.move_to((self.margins + (i * self.tick_spacing)), (self.height - 75))
            self.context.line_to((self.margins + (i * self.tick_spacing)), (self.height - 100))
            self.context.set_source_rgba (0, 0, 0)
            self.context.stroke()

    def add_xaxis(self):
        '''Add x-axis line.'''
        self.context.set_line_width(3)
        self.context.move_to((self.margins - 2), (self.height - 75))
        self.context.line_to((self.width - self.margins + 2), (self.height - 75))
        self.context.set_source_rgba (0, 0, 0)
        self.context.stroke()

    def add_ticklabels(self):
        '''Adds numerical sequence lengths to tick marks from add_ticks().'''
        for i in range(len(self.ticklabels)):
            self.context.set_source_rgb(0, 0, 0)
            self.context.set_font_size(20)
            self.context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            self.context.move_to((45 + (i * self.tick_spacing) - self.tickoffset[i]), (self.height - 50))
            self.context.show_text(self.ticklabels[i])
            self.context.stroke()

    def add_legend(self):
        '''Adds a legend of motif sequences and their respective colors in the sequences.'''
        # Creates box with border for legend contents
        legend_dim = ((self.width - 450), 10, (435), (len(self.motifs) * 25) + 15)
        self.context.rectangle(*legend_dim)
        self.context.set_source_rgba(0.8, 0.8, 0.8)
        self.context.fill()
        self.context.rectangle(*legend_dim)
        self.context.set_line_width(5)
        self.context.set_source_rgba(0, 0, 0)
        self.context.stroke()

        # Creating exon example and text
        self.context.rectangle((self.width - 430), 20, 40, 20)
        self.context.fill()
        self.context.set_line_width(5)
        self.context.move_to((self.width - 440), 40)
        self.context.line_to((self.width - 380), 40)
        self.context.stroke()
        self.context.set_font_size(30)
        self.context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.context.move_to((self.width - 360), 42)
        self.context.show_text("Exons")

        # Creating intron example and text with formating from exon entry
        self.context.move_to((self.width - 440), 70)
        self.context.line_to((self.width - 380), 70)
        self.context.stroke()
        self.context.move_to((self.width - 360), 80)
        self.context.show_text("Introns")

        # Adding motif colors and text
        for i, motif in enumerate(self.motifs):
            self.context.set_source_rgb(0, 0, 0)
            self.context.set_font_size(25)
            self.context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
            self.context.move_to((self.width - 180), (40 + (25 * i)))
            self.context.show_text(motif.seq)

            self.context.rectangle((self.width - 240), (20 + (25 * i)), (45), (20))
            self.context.set_source_rgba(*self.color_pairs[motif.seq])
            self.context.fill()

    def add_title(self):
        '''Add the input FASTA file name to the top of the output image.'''
        self.context.set_source_rgb(0, 0, 0)
        self.context.set_font_size(100)
        self.context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.context.move_to(25, 120)
        self.context.show_text(self.name)

    def add_axistitle(self):
        '''Add the axis title to the x axis for base position.'''
        self.context.set_source_rgb(0, 0, 0)
        self.context.set_font_size(25)
        self.context.select_font_face("Arial", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
        self.context.move_to(((self.width / 2) - 110), (self.height - 20))
        self.context.show_text("Number of Bases")
        self.context.stroke()

    def save_image(self):
        '''Saves the image as a PNG. File name is the name of the fasta file that was read in with all the sequences.'''
        self.surface.write_to_png(f"{self.name}.png")
        print(f"{self.name}.png saved!")

    def execute_all(self):
        '''Executes all methods for MotifsImage to output complete png.'''
        self.create_image()
        self.add_exons()
        self.add_readlabels()
        self.add_ticks()
        self.add_ticklabels()
        self.add_seqlines()
        self.add_motifs()
        self.add_legend()
        self.add_title()
        self.add_axistitle()
        self.add_xaxis()
        self.save_image()


##################################################################################################################################################################################################################################
# Creating Image
##################################################################################################################################################################################################################################
'''
The following lines of code will create an image from files input on the command line via argparse.
'''

if __name__ == "__main__":

    args = get_args()

    # Lists that will hold transcript and motif objects from input files
    transcripts_list = parse_fasta(args.fasta)
    motifs_list = parse_motif_txt(args.motifs)

    # Executing transcript method to find motif positions for all motifs by transcript
    for tran in transcripts_list:
        tran.find_motifs(motifs_list)

    # Capturing input FASTA file name for naming mm_image object
    file_name = re.findall(r"[^\/\\]+(?=\.fasta|\.fa)", args.fasta)[0]

    # Creating image object and executing all methods to output png image
    mm_image = MotifsImage(file_name, transcripts_list, motifs_list)
    mm_image.execute_all()
