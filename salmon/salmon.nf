#!/usr/bin/env nextflow

/*
  Copyright (c) 2021, icgc-argo-rna-wg

  Permission is hereby granted, free of charge, to any person obtaining a copy
  of this software and associated documentation files (the "Software"), to deal
  in the Software without restriction, including without limitation the rights
  to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
  copies of the Software, and to permit persons to whom the Software is
  furnished to do so, subject to the following conditions:

  The above copyright notice and this permission notice shall be included in all
  copies or substantial portions of the Software.

  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
  IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
  OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
  SOFTWARE.

  Authors:
    DabinJeong
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.1.0'  // package version

container = [
    'ghcr.io': 'icgc_test'
]
default_container_registry = 'ghcr.io'
/********************************************************************/


// universal params go here
params.container_registry = ""
params.container_version = "latest"
params.container = ""

params.cpus = 5
params.mem = 8 // GB
params.publish_dir = ""  // set to empty string will disable publishDir


// tool specific parmas go here, add / change as needed
params.input_file = "${baseDir}/tests/input/*_{1,2}.fastq.gz"
params.referenceSeq = "${baseDir}/input/gencode.v37.transcripts.fa"
params.annotation = "${baseDir}/input/*.gtf"
params.outdir = "${baseDir}/tests/expected/"
//params.output_pattern = "*"  // output file name pattern


inp_ch = Channel.fromFilePairs(params.input_file).ifEmpty{exit 1,"Fastq sequence not found: ${params.input_file}"}

process salmon {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:  // input, make update as needed
    tuple val(id), path(reads)
    path refSeq 
    path annotation 

  output:  // output, make update as needed
    path "out_dir/${id}.salmon", emit: output_file

  script:
    // add and initialize variables here as needed

    """
    python3 ${baseDir}/salmon.py \
      --referenceSeq ${refSeq} \
      --index salmon_index \
      --threads ${params.cpus} \
      --read1 ${reads[0]} \
      --read2 ${reads[1]} \
      --output-dir out_dir \
      --mem ${params.mem} > ${id}_expression_counting.log 2>&1
   
   awk 'FS=OFS="\t" {if (\$1!~/#/) print \$9}' $annotation |grep "gene_id"|grep "transcript_id"|awk -F" |; " -vOFS="\t" '{print \$4, \$2}'\
      > "${annotation}.tx2gene"
    
   Rscript ${baseDir}/tx2gene.R --input out_dir/quant.sf --tx2gene "${annotation}.tx2gene" --tool salmon --output "out_dir/${id}.salmon" 
   """
}

process norm{
  input:
    file readCounts

  output:
    publishDir "${params.outdir}"
    file "${readCounts}.norm"

  script:

    """
    python3 ${baseDir}/normalize.py $readCounts "${readCounts}.norm"
    """
}

// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  salmon(
    inp_ch,
    params.referenceSeq,
    params.annotation
  )
  norm(
   salmon.out
  )
}