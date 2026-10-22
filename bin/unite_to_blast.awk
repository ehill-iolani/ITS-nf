#!/usr/bin/awk -f
# UNITE general FASTA release -> a BLAST-ready fasta plus a lineage table.
#
# UNITE's headers are
#   >Fusarium_oxysporum|MH865823|SH1166021.10FU|reps|k__Fungi;p__Ascomycota;...;g__Fusarium;s__Fusarium_oxysporum
# which makeblastdb -parse_seqids rejects ('|' is its id separator), and whose
# lineage would otherwise end up in BLAST's stitle, where build_report.py, pcoa.py
# and the frontend all read the first two words as the species. So each record
# becomes
#   >MH865823 Fusarium oxysporum SH1166021.10FU
# and its lineage goes to a separate table (-v taxmap=<path>), one line per record:
#   MH865823 <tab> kingdom <tab> phylum <tab> class <tab> order <tab> family <tab> genus
# which build_report.py reads (--unite-taxmap) to write the abundance table's
# rank columns.
#
# The label (the first two words) follows 16S-nf's SILVA rule, which the
# frontend already understands: a named species is "Genus species"; a record
# UNITE couldn't name to species is "Genus sp.", or when the genus is unknown
# too, its lowest named rank with that rank in brackets ("Agaricaceae (family)").
# UNITE fills unknown ranks with "<higher taxon>_<rank>_Incertae_sedis" or
# "unidentified"; those count as no name.
#
# Plain awk (mawk) only: this runs in the stock ubuntu image.

BEGIN {
    FS = "|"
    split("kingdom phylum class order family genus", RANKS, " ")
    PREFIX["k"] = "kingdom"; PREFIX["p"] = "phylum"; PREFIX["c"] = "class"
    PREFIX["o"] = "order"; PREFIX["f"] = "family"; PREFIX["g"] = "genus"
    if (taxmap == "") {
        print "unite_to_blast.awk: pass -v taxmap=<path>" > "/dev/stderr"
        exit 1
    }
}

function placeholder(name) {
    return name == "" || name ~ /Incertae_sedis$/ || name ~ /^unidentified/
}

/^>/ {
    if (NF < 5) {
        print "unite_to_blast.awk: not a UNITE header: " $0 > "/dev/stderr"
        exit 1
    }
    acc = $2
    sh = $3
    split("", lin)
    species = ""
    n = split($5, levels, ";")
    for (i = 1; i <= n; i++) {
        p = substr(levels[i], 1, 1)
        name = substr(levels[i], 4)
        if (p == "s") species = name
        else if (p in PREFIX && !placeholder(name)) lin[PREFIX[p]] = name
    }

    if (species != "" && species !~ /_sp$/ && !placeholder(species)) {
        label = species
        gsub(/_/, " ", label)
    } else if ("genus" in lin) {
        label = lin["genus"] " sp."
    } else {
        label = "Unclassified"
        for (i = 5; i >= 1; i--) {
            if (RANKS[i] in lin) { label = lin[RANKS[i]] " (" RANKS[i] ")"; break }
        }
    }

    print ">" acc " " label " " sh
    row = acc
    for (i = 1; i <= 6; i++) row = row "\t" ((RANKS[i] in lin) ? lin[RANKS[i]] : "")
    print row > taxmap
    next
}

{ print toupper($0) }
