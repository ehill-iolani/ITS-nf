process BLAST_TAX {
    tag "$sample"
    label 'process_low'
    container "${params.container_registry}/biocontainers/blast:2.15.0--pl5321h6f7f691_1"
    publishDir(path: { "${params.outdir}/${sample}/07_taxonomy" }, mode: 'copy')

    input:
    tuple val(sample), path(consensus_fastas)
    path db_files
    val db_name

    output:
    tuple val(sample), path("${sample}.hits.tsv"), emit: hits

    script:
    // one blastn over all of a sample's consensus sequences, not one task per
    // cluster -- each task would otherwise stage the whole BLAST db (hundreds
    // of MB for SILVA) just to search a single sequence. File is sample-
    // prefixed since these get collected together in BUILD_REPORT.
    //
    // `awk 1` rather than cat so a fasta missing its trailing newline can't
    // run its last line into the next file's header
    """
    awk 1 ${consensus_fastas} > queries.fasta

    blastn -query queries.fasta -db ${db_name} \\
        -num_threads ${task.cpus} \\
        -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs stitle" \\
        -max_target_seqs 50 -evalue 1e-10 \\
        -out blast.tsv

    # 50 targets, not 16S-nf's 5: BUILD_REPORT ranks a cluster's hits by
    # identity and looks for near-ties among them, and BLAST's own cut to the
    # first N targets goes by score, which favours longer references -- with 5,
    # the closest UNITE records were often cut before BUILD_REPORT saw them.
    #
    # qcovs (how much of the consensus the hit covers, all HSPs together) is
    # what BUILD_REPORT's --min-qcov checks: ITS1 and ITS2 sit between
    # conserved 18S / 5.8S / 28S stretches, so an organism missing from the
    # reference can still hit something at high identity over the 5.8S alone.
    # It goes before stitle, which stays the last column.
    #
    # unmatched clusters (no hit) are flagged, not dropped -- a community
    # sample can easily contain strains that are simply absent from whatever
    # reference set was supplied, and that's a real (if uncertain) result,
    # not a pipeline failure
    #
    # qseqid here must match the id blastn used (the consensus fasta's header
    # token, "sample_clusterid") -- otherwise this row can't be joined back to
    # its sample in BUILD_REPORT. The FILENAME test (not NR == FNR) keeps this
    # right when blast.tsv is empty, i.e. when nothing at all hit.
    awk -F'\\t' '
        FILENAME == ARGV[1] { seen[\$1] = 1; next }
        /^>/ {
            id = substr(\$0, 2)
            sub(/[ \t].*/, "", id)
            if (!(id in seen)) printf "%s\\tNO_HIT\\tNA\\tNA\\tNA\\tNA\\tNA\\tflag_for_manual_review\\n", id
        }
    ' blast.tsv queries.fasta > no_hit.tsv

    cat blast.tsv no_hit.tsv > ${sample}.hits.tsv
    """
}
