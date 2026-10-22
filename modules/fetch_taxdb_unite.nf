process FETCH_TAXDB_UNITE {
    label 'process_low'
    // needs /bin/bash to exist (the Nextflow docker executor execs it
    // directly, regardless of the `shell` directive) -- ubuntu has it plus
    // tar/gunzip/awk already installed, only curl needs adding
    container "${params.dockerhub_registry}/library/ubuntu:22.04"
    storeDir params.tax_db_cache

    output:
    path "unite_10_euk.fasta", emit: fasta
    // each record's kingdom..genus, which BUILD_REPORT writes as the
    // abundance table's rank columns (see bin/unite_to_blast.awk)
    path "unite_10_euk_taxmap.tsv", emit: taxmap

    script:
    // UNITE 10.0 (2024-04-04) general FASTA release for *all eukaryotes*, not
    // the fungi-only one: ITS1/ITS4 also amplify plant (and some animal and
    // protist) ITS, and against a fungi-only reference those reads would come
    // back as poor fungal matches instead of as what they are. Dynamic
    // species-hypothesis representatives, singletons included. The same file
    // nf-core/ampliseq uses (https://doi.org/10.15156/BIO/2959334). The archive
    // also holds a "_dev" variant (untrimmed flanks), which is skipped.
    """
    apt-get update -qq && apt-get install -y -qq --no-install-recommends curl ca-certificates

    curl -fsSL --retry 3 -o unite.tgz \\
        "https://s3.hpc.ut.ee/plutof-public/original/1dda2021-4893-4f2f-b50e-87bfea795267.tgz"
    tar -xzf unite.tgz
    src=\$(ls sh_general_release_*.fasta | grep -v '_dev' | head -n 1)
    if [ -z "\$src" ]; then
        echo "ERROR: no sh_general_release_*.fasta in the UNITE archive" >&2
        exit 1
    fi

    unite_to_blast.awk -v taxmap=unite_10_euk_taxmap.tsv "\$src" > unite_10_euk.fasta
    rm -f unite.tgz sh_general_release_*.fasta

    if [ ! -s unite_10_euk.fasta ] || [ ! -s unite_10_euk_taxmap.tsv ]; then
        echo "ERROR: converted UNITE reference is empty" >&2
        exit 1
    fi
    """
}
