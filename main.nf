#!/usr/bin/env nextflow
/*
 * ITS-nf: ONT fungal ITS amplicon (ITS1/ITS4) classification for community
 * profiling. Forked from 16S-nf v0.3.1 -- the same clustering + consensus +
 * BLAST-taxonomy pipeline, with ITS references (UNITE, RefSeq ITS), an ITS
 * read-length window and a query-coverage check on BLAST hits. 16S-nf itself
 * has the same architecture as
 * edna-ont-nf (vsearch + spoa/racon/medaka + BLAST), loosely following
 * wf-metagenomics' stages (ingest -> filter -> classify -> abundance +
 * diversity report) but classifying via alignment against a targeted
 * reference (wf-metagenomics' "minimap2 mode", here via BLAST) rather than
 * Kraken2 -- a better fit once noisy per-read ONT data has been
 * consensus-polished, and it keeps the reference a plain FASTA a
 * non-specialist can swap out, same as edna-ont-nf's --taxdb.
 *
 * Emits the exact same abundance_table.tsv schema as edna-ont-nf, plus an
 * optional pcoa_coordinates.tsv for beta-diversity -- both consumed as-is by
 * the iolani-bioinformatics frontend's existing chart components.
 */

nextflow.enable.dsl = 2

include { WF_ITS                 } from './workflows/its.nf'
include { MERGE_FASTQ            } from './modules/merge_fastq.nf'
include { FETCH_TAXDB_UNITE      } from './modules/fetch_taxdb_unite.nf'
include { FETCH_TAXDB_REFSEQ_ITS } from './modules/fetch_taxdb_refseq_its.nf'

// ---- top-level params (override via -params-file or --flag) ----
// filtering/clustering/taxonomy defaults live in nextflow.config, not here --
// anything read directly inside an included module or workflow script
// (rather than only inside this file's own `workflow` block) must be set
// there to be reliably visible by the time that module/workflow script runs
params.input         = null   // path to samplesheet.csv (sample,fastq); fastq globs resolve relative to the launch dir, not the CSV's location
params.outdir        = "results"
params.fwd_primer    = null
params.rev_primer    = null
params.tax_db_source = 'unite'   // 'unite', 'refseq_its', or 'custom' (use --taxdb as-is) -- see helpMessage below
params.taxdb         = null   // path to an ITS reference sequences FASTA; required only when --tax_db_source custom. A BLAST db is built from this at runtime
params.metadata      = null   // optional CSV: sample,<arbitrary metadata columns> -- joined onto the PCoA output for coloring
params.help          = false

def helpMessage() {
    log.info """
    ITS-nf: ONT fungal ITS amplicon community classification
    ---------------------------------------------------------
    Usage:
      nextflow run main.nf --input samplesheet.csv --outdir results

    Required:
      --input       CSV: sample,fastq_path
      --taxdb       ITS reference sequences FASTA (required only when --tax_db_source custom)

    Reference database (--tax_db_source, default '${params.tax_db_source}'):
      unite        auto-download UNITE 10.0 general release, all eukaryotes, with lineage (cached under --tax_db_cache)
      refseq_its   auto-download NCBI RefSeq fungal ITS targeted loci (cached under --tax_db_cache)
      custom       use --taxdb as supplied (plain ">id Genus species ..." headers)

    Key optional:
      --fwd_primer / --rev_primer   primer sequences for cutadapt trimming (off by default)
      --min_len / --max_len / --min_qual   chopper filtering thresholds (defaults sized for ITS1-ITS4 amplicons)
      --enable_read_stats   read length/Q-score summary before vs. after filtering (default ${params.enable_read_stats})
      --cluster_id   vsearch identity threshold for clustering (default ${params.cluster_id})
      --min_cluster  minimum reads for a read-level cluster to get a consensus -- a compute/quality floor (default ${params.min_cluster})
      --spoa_max_reads  reads per cluster used for the spoa draft, 0 = all; racon still polishes with every read (default ${params.spoa_max_reads})
      --merge_id     identity at which per-cluster consensus sequences are merged, 0 to skip (default ${params.merge_id})
      --enable_medaka   use medaka-polished consensus instead of racon consensus (default ${params.enable_medaka})
      --blast_chunk_size  consensus sequences per BLAST task, 0 = one task per sample; same hits either way (default ${params.blast_chunk_size})
      --min_pident   BLAST %identity below which a hit is low-confidence (default ${params.min_pident})
      --min_qcov     % of the consensus a hit must cover, below which it is low-confidence (default ${params.min_qcov})
      --tie_pident   hits within this many % points of the best identity count as tied -> reported at genus level (default ${params.tie_pident})
      --min_abundance      merged clusters with fewer reads get low_abundance=true in the abundance table, not dropped (default ${params.min_abundance})
      --min_rel_abundance  ...or with a smaller fraction of the sample's clustered reads, 0 = off (default ${params.min_rel_abundance})
      --metadata     CSV: sample,<metadata columns> -- joined onto the PCoA coordinates for coloring by group/site/etc.
    """.stripIndent()
}

def samplesheetToChannel(path) {
    Channel
        .fromPath(path, checkIfExists: true)
        .splitCsv(header: true)
        .map { row ->
            // checkIfExists doesn't catch a glob that matches zero files, so check explicitly --
            // otherwise this silently produces an empty file list downstream instead of failing
            def fq = file(row.fastq, checkIfExists: true)
            def files = fq instanceof List ? fq : [fq]
            if (files.isEmpty()) {
                error "Samplesheet row '${row.sample}': no files matched '${row.fastq}' (glob paths resolve relative to the launch directory, not the samplesheet's location)"
            }
            tuple(row.sample, fq)
        }
}

workflow {
    if (params.help) {
        helpMessage()
        exit 0
    }

    if (!(params.tax_db_source in ['unite', 'refseq_its', 'custom'])) {
        log.error "Unknown --tax_db_source '${params.tax_db_source}' -- expected 'unite', 'refseq_its', or 'custom'"
        exit 1
    }

    // a missing input is a failed run, not a help request: exiting 0 here (as
    // 16S-nf v0.3.1 still does) reads as success to whatever launched the run
    if (!params.input || (params.tax_db_source == 'custom' && !params.taxdb)) {
        helpMessage()
        log.error !params.input ? "--input is required" : "--taxdb is required with --tax_db_source custom"
        exit 1
    }

    // --taxdb is only consulted for 'custom'; the two named sources fetch
    // (and cache under --tax_db_cache) their own reference fasta instead.
    // Only UNITE has a lineage table; the other sources pass a placeholder for
    // the same fixed-arity reason as NO_METADATA below
    taxmap_ch = Channel.fromPath("${projectDir}/assets/NO_TAXMAP")
    if (params.tax_db_source == 'unite') {
        FETCH_TAXDB_UNITE()
        taxdb_fasta_ch = FETCH_TAXDB_UNITE.out.fasta
        taxmap_ch = FETCH_TAXDB_UNITE.out.taxmap
    } else if (params.tax_db_source == 'refseq_its') {
        FETCH_TAXDB_REFSEQ_ITS()
        taxdb_fasta_ch = FETCH_TAXDB_REFSEQ_ITS.out.fasta
    } else {
        taxdb_fasta_ch = Channel.fromPath(params.taxdb, checkIfExists: true)
    }

    // PCOA's metadata input is optional -- when unset, pass a placeholder
    // file instead of an empty channel so the process (which has a fixed
    // two-input signature) always receives something to stage; the script
    // checks the staged filename to tell the two cases apart.
    metadata_ch = params.metadata
        ? Channel.fromPath(params.metadata, checkIfExists: true)
        : Channel.fromPath("${projectDir}/assets/NO_METADATA")

    WF_ITS(samplesheetToChannel(params.input), taxdb_fasta_ch, taxmap_ch, metadata_ch)
}

// dev/debug entry point -- run a single step in isolation, e.g.:
//   nextflow run main.nf -entry MERGE_ONLY --input samplesheet.csv -profile docker
workflow MERGE_ONLY {
    if (!params.input) {
        log.error "MERGE_ONLY requires --input samplesheet.csv"
        exit 1
    }
    MERGE_FASTQ(samplesheetToChannel(params.input))
}
