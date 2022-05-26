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
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/rna-seq-expression-counting.salmon'
]
default_container_registry = 'ghcr.io'
/********************************************************************/


// universal params go here
params.container_registry = ""
params.container_version = ""
params.container = ""

params.cpus = "5"
params.mem = "8" // GB
params.publish_dir = ""  // set to empty string will disable publishDir


// tool specific parmas go here, add / change as needed
//params.input_file = "NO_FILE_1/NO_FILE_1_{1,2}"
//params.referenceSeq = "NO_FILE_"
//params.annotation = "NO_FILE_3"
//params.outdir = ""
params.input_file="tests/input/sample_01_{1,2}.test.fastq.gz"
params.referenceSeq = "tests/input/gencode.v37.transcripts.fa"
params.annotation="tests/input/test.gtf"
params.output_pattern = "sample_01.test"  // output file name pattern


inp_ch = Channel.fromFilePairs(params.input_file).ifEmpty{exit 1,"Fastq sequence not found: ${params.input_file}"}

process salmon {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:  
    tuple val(id),path(reads)
    path referenceSeq
    path annotation 

  output: 
    publishDir
    tuple file("${params.output_pattern}.gene.salmon"),file("${params.output_pattern}.transcripts.salmon")

  script:
    // add and initialize variables here as needed

    """
    python3 /tools/salmon.py \
      --referenceSeq ${referenceSeq} \
      --index salmon_index \
      --read1 ${reads[0]} \
      --read2 ${reads[1]} \
      --output-dir out_dir \
      --threads ${params.cpus} \
      --mem ${params.mem} > ${params.output_pattern}_expression_counting.log 2>&1
   
   awk 'FS=OFS="\t" {if (\$1!~/#/) print \$9}' ${annotation} |grep "gene_id"|grep "transcript_id"|awk -F" |; " -vOFS="\t" '{print \$4, \$2}'\
      > "${annotation}.tx2gene"
   
   cp out_dir/quant.sf ./"${params.output_pattern}.transcripts.salmon"    
  
   Rscript /tools/tx2gene.R --input out_dir/quant.sf --tx2gene "${annotation}.tx2gene" --tool salmon --output "${params.output_pattern}.gene.salmon" 
   """
}

// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  salmon(inp_ch, file(params.referenceSeq), file(params.annotation))
  //salmon(inp_ch, params.annotation)
}
