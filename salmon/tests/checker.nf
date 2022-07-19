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
    'ghcr.io': 'ghcr.io/icgc-argo-workflows/rna-seq-expression-counting.salmon'
]
default_container_registry = 'ghcr.io'
/********************************************************************/

// universal params
params.container_registry = ""
params.container_version = ""
params.container = ""

// tool specific parmas go here, add / change as needed
params.input_file = ""
//params.index = "salmon_index"
params.referenceSeq = ""
params.annotation_file = ""
params.expected_output = ""

include {salmon} from '../salmon' 


process file_smart_diff {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"

  input:
    tuple path(output_file1),path(output_file2)
    path expected_file1
    path expected_file2 

  output:
    stdout()

  script:
    """
    # Note: this is only for demo purpose, please write your own 'diff' according to your own needs.
    # in this example, we need to remove date field before comparison eg, <div id="header_filename">Tue 19 Jan 2021<br/>test_rg_3.bam</div>
    # sed -e 's#"header_filename">.*<br/>test_rg_3.bam#"header_filename"><br/>test_rg_3.bam</div>#'
    for f in ${output_file1} ${expected_file1} 
        do     
        awk 'NR==1 {print};NR>1 {printf("%s %d %d\\n", \$1, int(\$2/10)*10, int(\$3/10)*10 )}' \$f > \$f".round" 
    done            

   for f in ${output_file2} ${expected_file2}
        do 
        awk 'NR==1 {print};NR>1 {printf("%s %d %d %d %d\\n", \$1, int(\$2/10)*10, int(\$3/10)*10, int(\$3/10)*10, int(\$3/10)*10 )}' \$f > \$f".round"
   done     
     
   diff ${output_file1}.round ${expected_file1}.round \
      && ( echo "Test PASSED" && exit 0 ) || ( echo "Test FAILED, output file mismatch." && exit 1 )
   
   diff ${output_file2}.round ${expected_file2}.round \
      && ( echo "Test PASSED" && exit 0 ) || ( echo "Test FAILED, output file mismatch." && exit 1 ) 
   """
}


workflow checker {
  take:
    inp_ch 
    annotation
    referenceSeq
    expected_output1
    expected_output2    

  main:
    salmon(
      inp_ch,
      referenceSeq,
      annotation 
    )
    file_smart_diff(
      salmon.out,
      expected_output1,
      expected_output2
    )
}


workflow {
  inp_ch = Channel.fromFilePairs(params.input_file).ifEmpty{exit 1,"Fastq sequence not found: ${params.input_file}"}
  checker(
    inp_ch, 
    file(params.annotation_file),
    file(params.referenceSeq),
    file(params.expected_output1),
    file(params.expected_output2)
  )
}
