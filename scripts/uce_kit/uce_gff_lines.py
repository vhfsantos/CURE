#!/usr/bin/env python
# coding: utf-8

# 24Oct2018 JBH
# given a uce data file in format shown below and a gff file, show gff lines matching uce pos for each uce
# also shows a summary line per uce including intergenic
# presumes uce data has mapped uces to a single point represented by its start pos:
#       uce-4323	NC_006088.4	2744945

import sys, re, time, os.path
from pprint import pprint

class Debug():
    debugging = False
    
    @staticmethod
    def print_scaf_counts(scaf_map, scaf_list):
        if not Debug.debugging:
            return
        for s in range(len(scaf_list)):
            scaf = scaf_list[s]
            sys.stderr.write("{} {}\n".format(scaf, len(scaf_map[scaf])))

    @staticmethod
    def pprint_scaf_uces(scaf_map, scaf):
        if not Debug.debugging:
            return
        for u in scaf_map[scaf]:
            pprint(vars(u))


def remove_version_suffix(name): # e.g. NC_006088.4 to NC_006088
    return re.sub("\.[0-9]*$","", name)


def to_int(digits_str, errval=0):  # do not throw a ValueError exception if input happens to not be a string of just digits
    return int(digits_str) if digits_str.isdigit() else errval


def get_subfld(fld, sep=";", which=0): # returns first semi-colon delimited field by default
    return fld.split(sep)[which]


def sort_scafmap_by_uce_pos(scaf_map):
    for scaf in scaf_map:
        scaf_map[scaf].sort(key=lambda uceinfo: uceinfo.pos)
        

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

def create_uce_scaf_map(ucefile, remove_version = True):
    class Uce_info(object): pass

    # create dict for each scaffold (name in fld 1) and store uce info
    #  (uce name, start pos, empty to start list of matching gff lines) for each uce in scaffold
    scaf_map = {}  # holds info about each scaf's uces
    scaf_list = [] # scaf names in order of occurrence so we can loop in same order later
    for ln in ucefile:
        ln = ln.rstrip("\n")
        uceflds = ln.split("\t")
        uce_nm = uceflds[0]
        uce_pos = int(uceflds[2])
        uce_scaf = uceflds[1]
        if remove_version:
            uce_scaf = remove_version_suffix(uce_scaf) # get rid of version info in scaffold name
        
        if not uce_scaf in scaf_map:
            scaf_map[uce_scaf] = []
            scaf_list.append(uce_scaf)
            
        # store info about the use in the relevant scaffold map
        uce_info = Uce_info()
        uce_info.uce = uce_nm;
        uce_info.pos = uce_pos
        uce_info.gff_lns = []
        
        scaf_map[uce_scaf].append(uce_info)
    
    # in case the ucefile was not in uce position sorted order, this will take care of that
    sort_scafmap_by_uce_pos(scaf_map)

    return scaf_map, scaf_list


def map_uces_to_gff_lines(uce_file, gff_file, begpos_ix = 3, endpos_ix = 4, exclude_list = ["region"], remove_version = True):
    ucefile = open(uce_file, "r")
    gff = open(gff_file, "r")
    
    uces_by_scaf_map, scaf_list = create_uce_scaf_map(ucefile, remove_version)
    Debug.print_scaf_counts(uces_by_scaf_map, scaf_list)

    for gff_ln in gff:
        # validate line
        if len(gff_ln) < 3 or gff_ln[0] == "#":  # ignore empty or comment lines
            continue

        cur_flds = gff_ln.split("\t")
        if len(cur_flds) <= endpos_ix: # line too short
            continue
        
        ln_type = cur_flds[2]
        if ln_type in exclude_list:
            continue
        
        begpos = to_int(cur_flds[begpos_ix])
        endpos = to_int(cur_flds[endpos_ix])
        if begpos < 1 or endpos < 1:
            continue
        
        if begpos > endpos: # swap them so we always have smaller position followed by larger
            begpos, endpos = endpos, begpos
        
        # line looks good, see if its extent has any of its scaffold's uce start positions in it
        scaf = cur_flds[0]
        if remove_version:
            scaf = remove_version_suffix(scaf)
            
        if not scaf in uces_by_scaf_map:
            continue
        
        # there are uces on this scaffold: see if this gff line is in one of the uce's scopes.
        # given certain assumptions about the order of the gff lines this could be done more
        # efficiently, however since isoforms make the order a little more complicated we will
        # do a search through all the specific scaffold's uce space for each line.
        
        for u in uces_by_scaf_map[scaf]:
            if begpos <= u.pos <= endpos: # this uce is the scope of this line
                u.gff_lns.append(gff_ln.rstrip("\n"))
                    
    return uces_by_scaf_map, scaf_list


def display_uce_gff_info(uces_by_scaf_map, scaf_list, show_summary = True, show_lines = False):
    show_summary = show_summary or not show_lines # if both are False we still want to show the summaries
    type_totals = {} # totals for exon_only, intron_only, exin_and_intron, or intergenic UCEs in the gff
    last_scaf = ""; last_pos = 0 # so we can output distance from last UCE as the 5th field
    
    num_uces = 0; num_gfflns = 0
    for scaf in scaf_list:
        for u in uces_by_scaf_map[scaf]:
            num_uces += 1; num_gfflns += len(u.gff_lns)
            distance = "" if scaf != last_scaf else u.pos-last_pos
            last_scaf = scaf; last_pos = u.pos
            uce_info = "{}\t{}\t{}\t".format(u.uce, scaf, u.pos)
            if len(u.gff_lns) > 0:
                if show_summary: # uce name, scaf, begpos, all the gff line types associated with this uce
                    typ = ""
                    for l in u.gff_lns:
                        flds = l.split("\t")
                        typ += flds[2]
                        if flds[2] == "gene":
                            gid = get_subfld(flds[8])
                            if gid != "":
                                typ += "(" + gid + ")"
                        typ += " "
                    sys.stdout.write("{}{}\t{}\t{}\n".format(uce_info, TypeUtil.shorthand_str(typ), distance, typ))
            
                if show_lines: # show uce name then gff line
                    for ln in u.gff_lns:
                        sys.stdout.write("{}\t{}\n".format(u.uce, ln))
            else:
                typ = "intergenic"
                sys.stdout.write("{}N\t{}\t{}\n".format(uce_info, distance, typ))
                
            TypeUtil.inc_totals_from_types(typ)
                                
    return num_uces, num_gfflns


def map_and_display(uce_filename, gff_filename, exclude_list, show_summary, show_lines):
    start = time.time()
    
    # gather up a map per scaffold of the uces in that scaffold and the gff lines overlapping each such uce
    uces_by_scaf_map, scaf_list = map_uces_to_gff_lines(uce_filename, gff_filename, exclude_list = exclude_list)
    
    # display the info based on defaults or user preferences that override them
    num_uces, num_gfflns = display_uce_gff_info(uces_by_scaf_map, scaf_list, show_summary, show_lines)

    # show what was done and how long it took    
    duration = (time.time() - start) / 1000 * 1000
    m, s = divmod(duration, 60); h, mrem = divmod(m, 60)
    
    sys.stderr.write("\r{} uces {} gff lines ({}m{:.3f}s)\n".format(num_uces, num_gfflns, int(m), s))
    sys.stderr.write("{}\n".format(TypeUtil.display_str()))


def usage(exit_code = 1):
    msg = """
    usage: uce_gff_lines.py <uce_name_pos_file> <gff_file> [-lines [-nosummary]] [<gff_line_type_to_exclude> ...]
    
    Input is a file with UCE info lines and a gff file, preferably with introns added
    (for this use you can use add_intron_to_gff.py or other programs).
    Each UCE info line should have 3 tabbed fields e.g: uce-4323	NC_006088.4	2744945
    where the 3rd field is the start position of the uce (end position is optional).
    
    Default output shows a summary for each UCE info line showing it, the number of gff lines
    and the type of each gff line overlapping the uce start position.
    
    If you use the -lines option, it also outputs each gff line that refers to the UCE
    outputting the UCE name prefixed as the first tabbed field of the line. When using
    the -lines option you can suppress summaries with -nosummary.
    
    Non-hyphen command line terms are considered types of gff lines to exclude from consideration.
    This might be CDS or mRNA. Comparisons are case sensitive. Term "region" excluded by default.
    
    The intergenic UCEs are shown in both cases. You can screen out any summary lines, including
    intergenic, and just retain the gff lines by piping output to: awk '! ($3~/^[0-9]+$/)' or you
    can remove gff_lines retaining only the summary by piping to: awk '($3~/^[0-9]+$/)'

"""
    sys.stderr.write(msg)
    sys.exit(exit_code)

    
def getoptions(argv, min_args, exit_code = 1):   
    if len(argv) < min_args:
        usage(exit_code)
        
    def checkfile(fname):
        if not os.path.isfile(fname):
            sys.stderr.write("File not found: {}\n".format(fname))
            sys.exit(exit_code)
    
    class options: pass
    options.uce_filename = argv[1]; checkfile(options.uce_filename)
    options.gff_filename = argv[2]; checkfile(options.gff_filename)
    options.show_summary = True
    options.show_lines = False
    options.excludes = ["region", "Region", "REGION"]
    
    
    # handle the options after the 2 file names. file names must be in those positions.
    for op in range(3,len(argv)):
        arg = argv[op]
        if arg[:3] == "-li":
            options.show_lines = True
        elif arg[:5] == "-nosu":
            options.show_summary = False
        elif arg[0] != '-':
            options.excludes.append(arg)
        elif arg == "-debug":
            Debug.debugging = True
        else:
            sys.stderr.write("invalid option: {}\n".format(arg))
            
    return options

    
def main(argv):  # pass in argv so we can load this file and call this as uce_gff_lines.main(sys.argv) from another python file                 
    ops = getoptions(argv, 3) # need at least 3 args: prog_name, uce_filename, gff_filename
    map_and_display(ops.uce_filename, ops.gff_filename, ops.excludes, ops.show_summary, ops.show_lines)


if __name__ == '__main__':
    from signal import signal, SIGPIPE, SIG_DFL
    signal(SIGPIPE,SIG_DFL) # don't complain if pipe closes output (head or less commands will cause this)
    
    try:
        main(sys.argv)
    except (KeyboardInterrupt, SystemExit): # don't do stack trace for Ctrl+C
        sys.exit(0)

# end of program
