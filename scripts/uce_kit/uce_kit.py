#!/usr/bin/env python
# coding: utf-8

# 26Oct2018 JBH
# routines for doing various things with uce files, including merging uce probes into single uces,
# filtering uces using a blat or blast m8 format file, also creating that file given uce file and
# genome file calling out to blat or blastn

# uce example record, probe 2 of 2:
"""
>uce-5_p2 |source:faircloth,probes-id:2683,probes-locus:5,probes-probe:2
AATTTCCTAGTTAAAACCCTCCCTTGCTGACAAGGGACTGAAAGAGTTTTAAATCACAGATGTAGAGTATCAAATGCAATAATGCTCTTGCAATAGTGCATTGAAGCCTCAATTAATTAA
"""

from __future__ import print_function # only using print for test function but make sure workins in python3 and python2
import sys, os.path
import re, operator, subprocess

# we parse the Uce fasta header line as a class method because we want to add all of the uce's probes
# to the same Uce object instance. So we can call parse_line and check to see if the uce
# name is already in a Uce object. That is how the UceList class uses this Uce list.
class Uce:
    overlap_bases = 60
    
    @classmethod
    def parse_line(cls, ln): # set class variables for the Uce line, e.g.: >uce-5_p2 |source:faircloth,probes-id:2683,probes-locus:5,probes-probe:2
        cls.init_vars()
        ln = ln.rstrip("\n")
        if len(ln) > 3:
            if ln[0] == ">":
                ln = ln[1:]
            flds = ln.split(" ")
            cls.uce_field = flds[0]
            cls.uce_comment = ln[(len(flds[0])+1) : ] # take rest of line after uce name, do this way so spaces in comment ok
            flds = cls.uce_field.split("_")
            cls.uce_name = flds[0]
            if len(flds) > 1:
                cls.uce_probe = flds[1]
            # print cls.uce_field; print cls.uce_name; print cls.uce_probe; print cls.uce_comment
            
        return cls.uce_name

    @classmethod
    def init_vars(cls):
        cls.uce_name = ""    # uce-5
        cls.uce_probe = ""   # p2
        cls.uce_comment = "" # |source:faircloth,probes-id:2683,probes-locus:5,probes-probe:2
        cls.uce_field = ""   # uce-5_p2
            
    # these are the instance methods
    def __init__(self, uce_name = None):  # call class method parse_line(ln) before creating object
        self.name = Uce.uce_name if uce_name == None else uce_name
        self.field = Uce.uce_field
        self.comment = Uce.uce_comment
        self.probes = [] # each entry is a list of probe name and probe sequence
        self.errormsg = ""
        
    def __len__(self):
        return len(self.probes)
    
    def __getitem__(self, ix_probe):  # get the probe tuple associated with ix_arg
        return self.probes[ix_probe]
    
    # define iter and next so can do for p in uce and get the probe lists
    def __iter__(self):
        self.ix_probe = 0
        return self
    def __next__(self):
        if self.ix_probe >= len(self.probes):
            raise StopIteration
        else:
            prb = self.probes[self.ix_probe]
            self.ix_probe += 1
            return prb
    next = __next__
        
    def get_probe(self, probe_number = 1): # return sequence for probe, first probe 1 by default
        num_probes = len(self.probes)
        if probe_number <= num_probes:
            return self.probes[probe_number-1][1]
        else:
            return "none"
    
    def get_probe_comment(self, probe_number = 1): # return the comment associated with the probe
        num_probes = len(self.probes)
        if probe_number <= num_probes:
            return self.probes[(probe_number-1)][2]
        else:
            return ""
        
    def get_probe_info(self, probe_number = 1): # return the whole 3 component list of the probe
        num_probes = len(self.probes)
        if probe_number <= num_probes:
            return self.probes[(probe_number-1)]
        else:
            return []
        
    def add_probe(self, sequence, probe = None, comment = None):
        probe = self.uce_probe if probe == None else probe
        comment = self.uce_comment if comment == None else comment
        self.probes.append([probe, sequence, comment])
        self.merge_probes() # merge them as they come in
        return self
        
    def merge_probes(self, overlap = 0): # probes overlap, usually 60 bases, merge them into a single sequence and return it
        overlap = self.overlap_bases if overlap == 0 else overlap
        self.errormsg = ""

        num_prbs = len(self.probes)
        self.sequence = "" if num_prbs < 1 else self.probes[0][1] # initialize with sequence of first probe
        for p in range(1, num_prbs):
            overlap_seq = self.sequence[-overlap:] # last 60 bases of the sequence
            probe_seq = self.probes[p][1]
            if probe_seq[:overlap] == overlap_seq: # probe looks good
                self.sequence += probe_seq[overlap:]
            else:
                self.errormsg = "{}: first {} bases of probe {} did not match last {} bases of probe {}.".format(self.uce_name, overlap, p+1, overlap, p)
                break 
        return self
    
    def write_uce_record(self, show_errs = True):
        num_probes = len(self.probes)
        sys.stdout.write( ">{} probes:{} {}\n{}\n".format(self.name, num_probes, self.get_probe_comment(), self.sequence) )
        if show_errs and self.errormsg != "":
            sys.stderr.write( "{} {} probes {}\n".format(self.name, num_probes, self.errormsg) )
            
    def write_probe_records(self):
        for p in self.probes:
            sys.stdout.write( ">{}_{} {}\n{}\n".format(self.name, p[0], p[2], p[1]))

# end Uce class definition

class UceList:
    def __init__(self):
        self.uce_map = {}
        self.uce_list = []
        
    def __len__(self):
        return len(self.uce_list)
    
    def __contains__(self, item):  # this is called for the "in" and "not in" membership tests
        return item in self.uce_map
    
    def __getitem__(self, arg):  # so can do container[2] or container["uce-12"] and get appropriate uce object returned
        name = self.uce_list[arg] if str(arg).isdigit() else arg
        return self.uce_map[name]

    def __iter__(self):  # implement __iter__ and __next__ to make this an iterable
        self.ix_uce_list = 0
        return self

    def __next__(self): # Python 3, Python 2 will call next():
        if self.ix_uce_list >= len(self.uce_list):
            raise StopIteration
        else:
            uce_nm = self.uce_list[ self.ix_uce_list ]
            self.ix_uce_list += 1
            return self.uce_map[uce_nm] # return object associated with the next uce_nm in the list
        
    next = __next__ # assign Python 2's next to Python 3's __next__ function, so they'll both work
    
    def add_probe_info(self, uceln, seq = None):
        uce_nm = Uce.parse_line(uceln) # use Uce classes to parse the line and set class variables
        if not uce_nm in self.uce_map and len(uce_nm) > 0: # first time we've seen this uce's probe, add it to our list and map
            self.uce_list.append(uce_nm)
            self.uce_map[uce_nm] = Uce(uce_nm) # create a Uce object and store it in the map
            
        if seq != None: # can do it in one go if we have the seq or call add_probe_sequence later when seq is known
            self.add_probe_sequence(uce_nm, seq)
        
        return uce_nm
    
    def add_probe_sequence(self, uce_nm, probe_seq):
        obj = None
        if uce_nm in self.uce_map:
            uce_obj = self.uce_map[uce_nm]
            obj = uce_obj.add_probe(probe_seq)
        return obj
    
    def load_file(self, uce_fname):
        if not checkfile(uce_fname):
            return None # None signals a problem
            
        ucefile = open(uce_fname, "r")
        sequence = "" # sequence of the current fasta record
        for ln in ucefile:
            ln = ln.rstrip("\n")
            if len(ln) < 3 or ln[0] == "#":
                continue
            if ln[0] == ">": # next fasta record started
                if sequence != "": # add sequence for the last fasta record
                    self.add_probe_sequence(uce_nm, sequence)
                    
                sequence = ""               
                uce_nm = self.add_probe_info(ln)
            else:
                sequence += ln
            
        if uce_nm != "" and sequence != "": # add the last probe sequence
            self.add_probe_sequence(uce_nm, sequence)
        return self
    
    def write_uce_records(self, show_errs = True):  # output merged probe sequences as fasta records
        for uce in self:
            uce.write_uce_record(show_errs)
    
#end UceList class definition


class TypeUtil:  # container for the routines that map uce type string types (exon, intron, etc) to other strings for display
    type_totals = {}
    
    @classmethod
    def inc_totals_from_types(cls, uce_typs): # e.g.: gene(ID=gene126) mRNA intron mRNA intron mRNA intron mRNA intron
        shorthand = cls.shorthand_str(uce_typs)
        cls.inc_totals(shorthand)

    @classmethod
    def display_str(cls): # type_totals map has shorthand str for key, "EI", "E", "I", "N", and count for value
        tss = cls.type_totals.copy()  # so we can delete items as we convert them
        cls.display = ""
        cls.total = 0
    
        def add_count(shrt, msg): # e.g., tss["E"] counts uces with only exon overlaps
            if shrt in tss:
                cls.total += tss[shrt]
                cls.display += msg.format(tss[shrt])
                del tss[shrt]
            
        add_count("E", "{} exon-only ")
        add_count("I", "{} intron-only ")
        add_count("EI","{} exon-and-intron ")
        add_count("N", "{} intergenic ")
            
        # now handle any that we didn't expect to occur
        for rare in tss:
            nm = rare if rare != "" else "unclassified"
            cls.display += "{} '{}' ".format(tss[rare], nm)
            cls.total += tss[rare]
            
        cls.display = "{} total: {}".format(cls.total, cls.display.rstrip(" "))
        return cls.display
        
    @classmethod
    def clear_totals(cls):
        type_totals.clear()
        
    # helper methods
    @classmethod
    def inc_totals(cls, shorthand): # "E", "I", "EI" is shorthand
        cls.type_totals[shorthand] = (cls.type_totals[shorthand] + 1) if shorthand in cls.type_totals else 1

    @staticmethod
    def shorthand_str(uce_typs): # we append a list to the totals_list with 'E', 'I', 'N' as needed
        uce_totals_counts = TypeUtil.uce_type_counts(uce_typs) # counts for the current uce
        typ_str = ""  # will have an E for exon, I for Intron and N for iNtergenic. so one with exon and intron would be "EI"
        if uce_totals_counts["exon"] > 0:
            typ_str += "E"
        if uce_totals_counts["intron"] > 0:
            typ_str += "I"
        if uce_totals_counts["intergenic"] > 0:
            typ_str += "N"
        return typ_str
    
    @staticmethod            
    def uce_type_counts(uce_typs):  # given the string with the uce gff line types, return a map of their counts
        uce_totals = {} # total for the current uce
        uce_totals["exon"]=0; uce_totals["intron"]=0; uce_totals["intergenic"]=0
        for w in uce_typs.split(" "):
            if w in ["exon", "intron", "intergenic"]:
                uce_totals[w] += 1       
        return  uce_totals
    
# end class TypeUtil


def show_summary_file_totals(summary_file):
    # summary file is a tsv file where lines with the 3rd field being all digits
    # describe the types of gff lines associated with a single UCE in the last field
    if not checkfile(summary_file):
        return

    num_uce_lines = 0
    fh = open(summary_file, "r")
    for tsv in fh:
        flds = tsv.rstrip("\n").split("\t")
        if len(flds) < 5 or not str(flds[2]).isdigit():
            continue
        num_uce_lines += 1
        TypeUtil.inc_totals_from_types( flds[-1] ) # last field has the gff types string
        
    if num_uce_lines > 0:
        sys.stdout.write("{}\n".format(TypeUtil.display_str()))
    else:
        sys.stderr.write("Looking for tabbed fields where 3rd field is an integer. None found.\n")


def filter_uce_match_m8_file(m8_filename, pct_match = 99.0, len_match = 120, exclude_prefixes = [], fout = sys.stdout):
    if not checkfile(m8_filename):
        return None # None signals a problem

    # for those lines where pct_match and len_match criteria are met use first probe of uce
    # in m8 file and write out uce name and first position in scaffold where the UCE starts
    # also exclude any lines where the scaffold name in fld 2 has a prefix in exclude_prefixes
    uce_info = [] # each entry will be a 3 value tuple
    uce_nm_map = {} # 10Jul2019 to keep track if we have already seen this uce
    last_nm = ""

    fh = open(m8_filename, "r")
    for tsv in fh:
        flds = tsv.rstrip("\n").split("\t")
        if len(flds) < 11:
            sys.stderr.write("not enough fields in line: {}".format(len(flds)))
            break
        
        if float(flds[2]) < pct_match or int(flds[3]) < len_match:
            continue

        probe_nm = flds[0] # consists of uce name and probe number e.g. uce-501_p1 or uce-36_6
        scaff_nm = flds[1]
        pos = flds[8] if flds[8] < flds[9] else flds[9]

        # JBH 10Jul2019 change method for determining uce to include (was based on _p1 which doesn't necessarily hold)
        ''' # old method
        uce_nm = re.sub("_p1$", "", probe_nm)
        if re.search("_p[0-9]+$", uce_nm): # it's not the first probe
            continue
        '''
        # new method, takes first one and ignores others with same prefix before _p[0-9]+
        uce_nm = re.sub("_p[0-9]+$", "", probe_nm)
        if uce_nm in uce_nm_map:
            continue
        uce_nm_map[ uce_nm ] = pos
        
        prefix_ok = True
        for p in exclude_prefixes:
            if re.search("^"+p, scaff_nm):
                prefix_ok = False
                break
        if not prefix_ok:
            continue
        
        # make sure we aren't trying to put in one that we just appended
        if last_nm == probe_nm:
            continue
        last_nm = probe_nm
        
        tup = (uce_nm, scaff_nm, int(pos), flds[2], flds[3])
        uce_info.append(tup)
     
    # sort by scaffold name in tup[1] then pos in tup[2]   
    uce_info.sort(key = operator.itemgetter(1, 2))
    for tup in uce_info:
        fout.write("{}\t{}\t{}\n".format(tup[0], tup[1], tup[2]))

def cmd_exists(cmd):
    return any(
        os.access(os.path.join(path, cmd), os.X_OK) 
        for path in os.environ["PATH"].split(os.pathsep)
    )


def check_pgm_and_ops(pgm_list, file1, file2):
    for pgm in pgm_list:
        if cmd_exists(pgm):
            break
    if not cmd_exists(pgm):
        sys.stderr.write("Program {} can't be found.\n".format(pgm_list[0]))
        return None
    
    checkfile(file1)
    checkfile(file2)   
    return pgm

def execute_cmdline(cmdline, outfile=None, msg_to_stderr=None):
    fh = None # this means use stdout    
    if outfile != None and outfile != "":
        fh = open(outfile, 'w')
    
    if msg_to_stderr != None:
        sys.stderr.write(msg_to_stderr)
        
    retcode = subprocess.call ( cmdline.rstrip().split(" "), stdout=fh )
    
    if fh != None:
        fh.close()
    return retcode
    

# call out to blat (or blatq) to do search of uce probes in genome fasta
# use blat setting minScore=100 to reduce cluttered hits fot the 120 base probes
def blat_uces(uces, genome, output):
    pgm = check_pgm_and_ops(["blat", "blatq"], uces, genome)
    if pgm == None:
        return -1
    
    if output == "": output = "stdout"
    settings = "-stepSize=5 -repMatch=100000 -out=blast8 -minScore=100"
    cmdline = "{} {} {} {} {}".format(pgm, settings, genome, uces, output)
    sys.stderr.write("Searching for uces: {}\n".format(cmdline))
    retcode = subprocess.call( cmdline.split(" ") )
    return retcode


# makeblastdb for the genome file in ops.filename if it does not exist.
# makeblastdb -db galGal6.fasta -dbtype nucl
# makes galGal6.fasta.nin galGal6.fasta.nhr & galGal6.fasta.nsq
# blastn -query genes.ffn -subject genome.fa -outfmt 6 -evalue 1e-50
def blast_uces(uces, genome, eVal, addtl_args = ""):
    pgm = check_pgm_and_ops(["blastn", "makeblastdb"], uces, genome)
    if pgm == None:
        return -1

    evalue = "-evalue " + eVal
    if addtl_args.find("-evalue ") > -1: # can't define same arg twice for blast utilities
        evalue = ""
    
    ext = [".nin",".nhr",".nsq"] # extensions of files created by makeblastdb added to the genome fasta filename
    db_exists = os.path.isfile(genome+ext[0]) and os.path.isfile(genome+ext[1]) and os.path.isfile(genome+ext[2])
    
    if not db_exists:
        cmdline = "makeblastdb -in {} -dbtype nucl".format(genome)
        sys.stderr.write(cmdline)
        retcode = subprocess.call ( cmdline.split(" "), stdout=sys.stderr )
        if retcode != 0:
            sys.stderr.write("error {} creating blast db\n".format(retcode))
            return
        sys.stderr.write("\n")
        
    cmdline = "blastn -query {} -db {} -outfmt 6 {} {}".format(uces, genome, evalue, addtl_args)
    sys.stderr.write("Searching for uces: {}\n".format(cmdline))
    retcode = subprocess.call ( cmdline.rstrip(" ").split(" ") )
    return retcode

def run_pipeline(ops, output_dir):
    def remove_ext(fname):
        lst = fname.split(".")
        fname = fname.split("/")
        fname = fname[len(fname)-1]
        if len(lst) < 2: return fname
        lst.pop()
        return ".".join(lst)
    def add_before_ext(fname, new_itm):
        lst = fname.split(".")
        if len(lst) < 2: return "{}.{}".format(fname, new_itm)
        lst.append(lst[-1]); lst[-2] = new_itm
        return ".".join(lst)
    
    uce = ops.filename; genome = ops.genome; gff = ops.gff
    checkfile(uce); checkfile(genome); checkfile(gff)
    
    # setting full path of this file to be used later...
    wd = os.path.abspath(__file__)
    wd = wd.split("uce_kit.py")[0]
        
    # Step 1 blat (or blastn) tetrapod_uces.fasta galGal4.fna >tetrapod_uces.galGal4_matches.m8
    m8_file = output_dir+"/{}.matches.m8".format(remove_ext(os.path.basename(uce)))
    sys.stderr.write("\nStep 1 of 4: ") # blat_uces() will write rest of line
    if cmd_exists("blat") or cmd_exists("blatq"):
        rc = blat_uces(uce, genome, m8_file)
    else: # could not find blat or blatq, run blast
        rc = blast_uces(uce, genome, ops.evalue, "-out " + m8_file)

    if rc != 0:
        return

    # Step 2 filter_tsv tetrapod_uces.galGal4_matches.m8 NT_ >galGal4_uce_locations.tsv
    uce_locs_name = "{}_uce_locations.tsv".format(output_dir+"/"+remove_ext(os.path.basename(genome)))
    filt_file = uce_locs_name
    addtl = "" if len(ops.filter_exclude_list) == 0 else " ".join(ops.filter_exclude_list) + " "
    sys.stderr.write("Step 2 of 4: filter_tsv {} {}>{}\n".format(m8_file, addtl, filt_file))
    fh_uce_locs = open(filt_file, 'w')
    filter_uce_match_m8_file(m8_filename = m8_file, exclude_prefixes = ops.filter_exclude_list, fout = fh_uce_locs)
    fh_uce_locs.close()

    # Step 3 add_introns_to_gff.py galGal4.gff >galGal4.with_introns.gff
    # next line modified by vhfsantos, 2021
    new_gff_name = output_dir+"/"+os.path.basename(gff) + ".with.introns"
    cmdline = wd+"add_introns_to_gff.py {}".format(ops.gff)
    retcode = execute_cmdline(cmdline, new_gff_name, "Step 3 of 4: {} >{}\n".format(cmdline, new_gff_name))
    
    # Step 4 uce_gff_lines.py galGal4_uce_locations.tsv galGal4_with_introns.gff >galGal4_uce_type_summary.txt
    # next line modified by vhfsantos, 2021
    summary_fname = "{}.uce_kit_summary".format(output_dir+"/"+os.path.basename(genome))
    cmdline = wd+"uce_gff_lines.py {} {}".format(uce_locs_name, new_gff_name)
    cmdline += " -lines" if ops.gff_lines else ""
    cmdline += " " + " ".join(ops.gff_exclude_list) if len(ops.gff_exclude_list) > 0 else ""
    retcode = execute_cmdline(cmdline, summary_fname, "Step 4 of 4: {} > {} \n\n".format(cmdline, summary_fname))
    
    sys.stderr.write("\n")
    return # from run_pipeline()

def test(): # testing code
    uce_container = UceList()
    uceln = ">uce-5_p1 |source:faircloth,probes-id:2682,probes-locus:5,probes-probe:1"
    uceseq = "AGCATCCTTAATATTTCCTTCCTTTCAGAAGCAAATAGAGCGTACCCTTATCTGAATGCTAATTTCCTAGTTAAAACCCTCCCTTGCTGACAAGGGACTGAAAGAGTTTTAAATCACAGA"

    uce_nm = uce_container.add_probe_info(uceln, uceseq) # adds uce if first probe for this one found

    uceln = ">uce-5_p2 |source:faircloth,probes-id:2683,probes-locus:5,probes-probe:2"
    uceseq = "AATTTCCTAGTTAAAACCCTCCCTTGCTGACAAGGGACTGAAAGAGTTTTAAATCACAGATGTAGAGTATCAAATGCAATAATGCTCTTGCAATAGTGCATTGAAGCCTCAATTAATTAA"

    uce_nm = uce_container.add_probe_info(uceln) # adds uce if first probe for this one found
    uce_obj = uce_container.add_probe_sequence(uce_nm, uceseq)
    
    print(uce_obj.sequence)
    print( "\nlen(uce_obj) shows {} probes in uce_obj {}".format(len(uce_obj), uce_obj.name) )
    print( "\nuce_obj[1] probe info: {}\n".format(uce_obj[1]) )

    print( "write probes with uce_obj.write_probe_records():" )
    uce_obj.write_probe_records()
    
    print( "\niterate over probes with 'for p in uce_obj:' and show probe name" )
    for p in uce_obj: print("{}_{}".format(uce_obj.name, p[0]))
    
    print("\nlen(uce_container):",len(uce_container))
    print("\nIndex by int uce_container[0]: \n", vars(uce_container[0]), "\n")
    print("Index by name uce_container[\"uce-5\"]\n", vars(uce_container["uce-5"]), "\n\n")
    print("'missing' in uce_container: ", 'missing' in uce_container)
    print("'uce-5' in uce_container: ", 'uce-5' in uce_container)
    print("\nuce_container['missing']", uce_container['missing'])  # this should blow up with KeyError and stack trace
    return

def pipeline_doc():
    doc = """
    Given the three files, (1) UCE file in fasta format, (2) genome fasta file, and (3) genome's gene annotations in gff format,
    this shows how to create a file that shows summaries, and optionally detail, of the UCE gene types. That is, whether a
    UCE is one or more of Exon, Intron or Intergenic type in the subject genome.
    
    Assume the three files above are named: tetrapod_uces.fasta, galGal4.fna, and galGal4.gff
    Here's an abbreviated pipeline, with $ referring to the command line prompt of the terminal:
    
    $ uce_kit.py blat tetrapod_uces.fasta galGal4.fna >tetrapod_uces.galGal4_matches.m8
    
    $ uce_kit.py filt tetrapod_uces.galGal4_matches.m8 NT_ >galGal4_uce_locations.tsv
    
    $ add_introns_to_gff.py galGal4.gff >galGal4.with_introns.gff
    
    $ uce_gff_lines.py galGal4_uce_locations.tsv galGal4_with_introns.gff >galGal4_uce_type_summary.txt
    
    You can also have the 4 steps run for you using:
    
    $ uce_kit.py run_pipeline tetrapod_uces.fasta galGal4.fna galGal4.gff -filt NT_
    
    The above gets the summary of the types with additional info such as gene ID.
    Here's the command line info for that last program, you can run it with no arguments for a description:

        uce_gff_lines.py <uce_name_pos_file> <gff_file> [-lines [-nosummary]] [<gff_line_type_to_exclude> ...]
        
    Also note that last argument to the uce_kit.py filt tool. It is NT_ and means to exclude all scaffolds with
    UCE matches to the genome that occur in scaffolds that have an NT_ prefix.
    
    Example summary output. Gff lines not included in this. Distance is from previous UCE on scaffold.
    Columns: UCE    Scaf            UCE pos         Type    Distance GFF type list:
    uce-6357        NC_006088       16361333        I       238419   gene(ID=gene299) mRNA intron mRNA intron mRNA intron 
    uce-6734        NC_006088       17315669        EI      954336   gene(ID=gene301) mRNA intron mRNA intron mRNA exon mRNA exon mRNA exon 
    uce-6265        NC_006088       17632542        E       316873   gene(ID=gene301) mRNA exon mRNA exon mRNA exon 
    uce-7988        NC_006088       17765381        N       132839   intergenic
    uce-4637        NC_006088       18152376        N       386995   intergenic
    uce-4779        NC_006088       18364070        N       211694   intergenic

 """       
    sys.stdout.write(doc)
    return


def checkfile(fname, exit_on_fail = True, exit_code = 1):  # call this for each argument expected to be an existing file
    if not os.path.isfile(fname):
        sys.stderr.write("File not found: {}\n".format(fname))
        if exit_on_fail:
            sys.exit(exit_code)
        return False
    return True

def usage(exit_code = 1):
    msg = """
    usage: uce_kit.py blat | blast | filter_tsv | summary_totals | merge | run_pipeline  [file1] <file2> <option1>... | pipeline_doc
    
    filt[er_tsv] <blast_or_blat_results> <prefix_to_exclude>... # filter the results to output the uce_info file needed by uce_gff_lines.py

    blat <uce_fasta_file> <genome_fasta_file>   # uses blat to find where uce's are in the genome file, outputs in tsv format (no indexing required)
    blast <uce_fasta_file> <genome_fasta_file>  # uses blastn to find where uce's are in the genome file, outputs in tsv format (db indexing done first)
    
    merge <uce_probe_fasta_file>                # merge individual uce probes into a single sequence and output as fasta records
    
    sum[mary_totals] <uce_gff_summary_file>     # one line total for the uce summary types from the output of uce_gff_lines.py
    
    pipe[line_doc]                              # short description of how to go from uce's and genomes to uce exon, intron, intergenic categorization
    
    run_pipeline <uces> <genome> <gff>          # given these three files, run the 4 step pipeline documented by uce_kit.py pipeline_doc. -filt excludes
     [-filt <ex1>...] [-excl <ex1>...] [-lines] #  items with prefix in filter Step 2, -excl exclude terms and -lines outputs gff lines in gff Step 4.

"""
    sys.stderr.write(msg)
    sys.exit(exit_code)

def getops(argv, min_args = 0, exit_code = 1):
    num_args = len(argv)-1 # argv[0] is program name, which we don't consider arg for this context
    if num_args < min_args:
        usage(exit_code)
        
    class options:
        action = argv[1]
    
    if options.action == "help":
        usage(exit_code)
    elif options.action == "test":
        return options
    elif options.action[:4] == "pipe":
        options.action = "pipeline_doc"
        return options
    
    if num_args < min_args+1: # need at least one more arg if this isn't calling for a test
        usage(exit_code)

    options.filename = argv[2]            
    options.evalue = "9e-40"
    if options.action == "blast" or options.action == "blat":
        options.output = ""   # empty str means write to stdout
        if num_args < 3:
            sys.stderr.write("\n    {} requires 2 files\n".format(options.action))
            usage(exit_code)
        options.file2 = argv[3]
        # options for blast
        options.addtl_args = ""
        for ix in range(4, len(argv)):
            options.addtl_args += argv[ix] + " "
            
    elif options.action[:4] == "filt":
        options.action = "filter_tsv"
        options.exclude_list = []
        for ix in range(3, len(argv)):
            arg = argv[ix]
            if arg[0] == '-' and (ix+1) < len(argv):
                parm = argv[ix+1]
                if arg[:4] == "-pct":
                    options.pct = float(parm)
                elif arg[:4] == "-len":
                    options.len = parm if parm.isdigit else "120"
            else:
                options.exclude_list.append(arg)
                
    elif options.action[:3] == "sum":
        options.action = "summary_totals"
        
    elif options.action[:3] == "run":
        options.action = "run_pipeline"
        if num_args < 4:
            sys.stderr.write("\n    {} requires 3 files\n".format(options.action))
            usage(exit_code)
        options.genome = argv[3]
        options.gff = argv[4]
        options.gff_exclude_list = []; options.filter_exclude_list = []
        excl_list = options.gff_exclude_list
        options.gff_lines = False
        output_dir = argv[5]
        
        for ix in range(6, len(argv)-1):
            arg = argv[ix]
            sys.stderr.write("{} :: {}\n". format(ix, arg))
            if arg[0] == '-': # -filt sets up for filter prefix exclusion list else the list is for the gff type to exclude
                if arg[:4] == "-lin": # add -lines to step 4 so gff lines included
                    options.gff_lines = True
                else:
                    excl_list = options.filter_exclude_list if arg[:5] == "-filt" else options.gff_exclude_list
            else:
                excl_list.append(arg)

                
    elif not options.action == "merge":
        sys.stderr.write("\n    \"{}\" is an invalid operation\n".format(options.action))
        usage(exit_code)

    return options, output_dir

def main(argv):  # pass in argv so we can load this file and call this from another python file                 
    ops, output_dir = getops(argv, min_args = 1, exit_code = 1) # exit_code is used if something wrong with the options
    
    if ops.action == "test":
        test()
    elif ops.action == "pipeline_doc":
        pipeline_doc()
    elif ops.action == "blat":
        blat_uces(ops.filename, ops.file2, ops.output)
    elif ops.action == "blast":
        blast_uces(ops.filename, ops.file2, ops.output, ops.evalue, ops.addtl_args)
    elif ops.action == "filter_tsv":
        filter_uce_match_m8_file(ops.filename, exclude_prefixes = ops.exclude_list)
    elif ops.action == "summary_totals":
        show_summary_file_totals(ops.filename)
    elif ops.action == "merge":
        uce_list = UceList()
        if uce_list.load_file(ops.filename):
            uce_list.write_uce_records()
    elif ops.action == "run_pipeline":
        run_pipeline(ops, output_dir)
    
    return

if __name__ == '__main__':
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL) # don't complain if pipe closes output (head or less commands will cause this)
    
    try:
        main(sys.argv)
    except (KeyboardInterrupt, SystemExit): # don't do stack trace for Ctrl+C
        sys.exit(0)

# end of program
