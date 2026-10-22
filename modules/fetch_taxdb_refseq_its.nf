process FETCH_TAXDB_REFSEQ_ITS {
    label 'process_low'
    // needs /bin/bash to exist (the Nextflow docker executor execs it
    // directly, regardless of the `shell` directive) -- ubuntu has it plus
    // gunzip already installed, only curl needs adding
    container "${params.dockerhub_registry}/library/ubuntu:22.04"
    storeDir params.tax_db_cache

    output:
    path "refseq_its.fasta", emit: fasta

    script:
    // NCBI's RefSeq Targeted Loci fungal ITS export: curated, mostly from type
    // material, fungi only and far smaller than UNITE. Already plain DNA fasta
    // with "accession Genus species ..." headers. It carries no lineage, so a
    // RefSeq run rolls up no higher than genus.
    //
    // Its titles end "ITS region; from TYPE material", and a ';' in a BLAST
    // title is how build_report.py, pcoa.py and the frontend recognise a SILVA
    // lineage -- so the headers' semicolons become commas, or every RefSeq hit
    // would be parsed as a lineage
    """
    apt-get update -qq && apt-get install -y -qq --no-install-recommends curl ca-certificates

    curl -fsSL --retry 3 -o fungi.fna.gz \\
        "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz"
    gunzip -c fungi.fna.gz | awk '/^>/ { gsub(/;/, ",") } { print }' > refseq_its.fasta
    rm -f fungi.fna.gz

    if [ ! -s refseq_its.fasta ]; then
        echo "ERROR: downloaded RefSeq ITS reference is empty" >&2
        exit 1
    fi
    """
}
