#!/usr/bin/env nextflow

/*
 * Pipeline parameters
 */


//Output directories


params.suppa_out="suppa.dir"


/*
 * Processes
 */


//suppa generate exon skipping events from model gtf
process suppa_ge {

        conda params.novelTranscripts_conda

        publishDir params.suppa_out, mode: 'symlink'

	input:
            path models_gtf

	output:
	    tuple path("${models_gtf}"), path("transcript_model_as_events_SE_strict.ioe")

	script:
	"""
        suppa.py generateEvents -i '$models_gtf' \
                                -o transcript_model_as_events \
                                -f ioe \
                                -e SE
	"""
}

//suppa calculate psi (relative abundance) for each transcript
process suppa_isoPsi {

        conda params.novelTranscripts_conda

        publishDir params.suppa_out, mode: 'symlink'

        input:
            tuple path(models_gtf), path(se_ioe)

        output:
            tuple path("${models_gtf}"), path("${se_ioe}"), path("transcript_model_isoform.psi")

        script:
        """
        suppa.py psiPerIsoform -g '$models_gtf' \
                               -e '${params.model_tpms}' \
                               -o transcript_model
        """
}

//convert se ioe file to bed file of cassette exons
process ioe2bed {

        conda params.novelTranscripts_conda

        publishDir 'files', mode: 'symlink'

        input:
            tuple path(models_gtf), path(se_ioe), path(iso_psi)

        output:
            tuple path("${models_gtf}"), path("cassetteExon.bed"), path("${iso_psi}")

        script:
        """
        Rscript '$projectDir'/scripts/ioe2bed.R -e '$se_ioe'
        """
}

//overlap cassette exon bed with gencode ref gtf exons using bedtools intersect -wao (null overlaps reported) 
process wao_intersect {

        conda params.novelTranscripts_conda

        publishDir 'files', mode: 'symlink'

        input:
            tuple path(models_gtf), path(ce_bed), path(iso_psi)

        output:
            tuple path("${models_gtf}"), path("ce_vs_ref.txt"), path("${iso_psi}")

        script:
        """
        bedtools intersect -wao -s -a '$ce_bed' -b <(gzcat '${params.refgtf}' | awk '\$3=="exon"') > ce_vs_ref.txt
        """
}

//putative cassette exon QC, associate with transcripts and write final output tables
process get_novel_transcripts {

        conda params.novelTranscripts_conda

        publishDir 'files', mode: 'symlink'

        input:
            tuple path(models_gtf), path(ce_vs_ref), path(iso_psi)

        output:
            tuple path("novelTranscripts.tsv"), path("novelCassetteExons.tsv")

        script:
        """
        Rscript '$projectDir'/scripts/getNovelTranscripts.R -c '$ce_vs_ref' \
                                                            -m '$models_gtf' \
                                                            -p '$iso_psi' \
                                                            -t '${params.model_tpms}'
        """
}

/*
 * Workflow
 */


workflow {

    //create input channel
    models_ch = channel.fromPath(params.models)
    
    suppa_ge(models_ch)

    suppa_isoPsi(suppa_ge.out)

    ioe2bed(suppa_isoPsi.out)

    wao_intersect(ioe2bed.out)

    get_novel_transcripts(wao_intersect.out)
}

