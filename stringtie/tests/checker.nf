#!/usr/bin/env nextflow

/*
  Copyright (c) 2022, icgc-argo-workflows

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

/*
 This is an auto-generated checker workflow to test the generated main template workflow, it's
 meant to illustrate how testing works. Please update to suit your own needs.
*/

/********************************************************************/
/* this block is auto-generated based on info from pkg.json where   */
/* changes can be made if needed, do NOT modify this block manually */
nextflow.enable.dsl = 2
version = '0.1.0'  // package version

container = [
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/expression-counting.stringtie'
]
default_container_registry = 'ghcr.io'
/********************************************************************/

// universal params
params.container_registry = ""
params.container_version = ""
params.container = ""

// tool specific parmas go here, add / change as needed
params.input_file = ""
params.annotation = ""
params.expected_output1 = ""
params.expected_output2 = ""

include { stringtie } from '../stringtie' 


process file_smart_diff {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"

  input:
    path output_file
    path expected_file1
    path expected_file2 

  output:
    stdout()

  script:
    """
    sed -e "s\/gene_id,.*$\/gene_id,readCounts\/g" ${output_file[0]} > ${output_file[0]}".diff" 
    sed -e "s\/gene_id,.*$\/gene_id,readCounts\/g" ${expected_file1}> ${expected_file1}".diff"
    
    sed -e "s\/transcript_id,.*$\/transcript_id,readCounts\/g" ${output_file[1]} > ${output_file[1]}".diff"    
    sed -e "s\/transcript_id,.*$\/transcript_id,readCounts\/g" ${expected_file2} > ${expected_file2}".diff"   
    
    diff ${output_file[0]}".diff" ${expected_file1}".diff" \
        && ( echo "Test PASSED" && exit 0 ) || ( echo "Test FAILED, output file of gene-level stringtie quantification mismatch." && exit 1 )
    
    diff ${output_file[1]}".diff" ${expected_file2}".diff" \
        && ( echo "Test PASSED" && exit 0 ) || ( echo "Test FAILED, output file of transcript-level stringtie quantification mismatch." && exit 1 )
    """
}


workflow checker {
  take:
    input_file
    annotation
    expected_output1
    expected_output2

  main:
    stringtie(
      input_file,
      annotation
    )

    file_smart_diff(
      stringtie.out,
      expected_output1,
      expected_output2
    )
}


workflow {
  checker(
    params.input_file,
    params.annotation,
    params.expected_output1,
    params.expected_output2
  )
}
