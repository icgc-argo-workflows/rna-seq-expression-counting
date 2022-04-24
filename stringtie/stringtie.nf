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


// universal params go here
params.container_registry = ""
params.container_version = ""
params.container = ""

params.cpus = 1
params.mem = 1  // GB
params.publish_dir = ""  // set to empty string will disable publishDir


// tool specific parmas go here, add / change as needed
params.input_file = "${baseDir}/tests/input/*.bam"
params.annotation = "${baseDir}/tests/input/*.gtf"
params.outdir = "${baseDir}/tests/expected/"
//params.output_pattern = "*"  // output file name pattern

inp_bam_ch = Channel.fromPath(params.input_file).map{ file->tuple(file.baseName, file) }.ifEmpty{exit 1, "bam file not found: ${params.input_file}"}


process stringtie {
  container "${params.container ?: container[params.container_registry ?: default_container_registry]}:${params.container_version ?: version}"
  publishDir "${params.publish_dir}/${task.process.replaceAll(':', '_')}", mode: "copy", enabled: params.publish_dir

  cpus params.cpus
  memory "${params.mem} GB"

  input:  // input, make update as needed
    tuple val(id), path(bam)
    path annotation

  output:  // output, make update as needed
    tuple val(id), file("${id}.stringtie"), file("${id}.stringtie.gtf")
    //path "out_dir/${params.output_pattern}", emit: output_file

  script:
    // add and initialize variables here as needed
    """
    mkdir -p out_dir

    python3 ${baseDir}/stringtie.py \
      -a ${annotation} \
      -bam ${bam} \
      -o ${id} 
    """
}

process stringtie_rc{

        input:
        tuple val(id), file(stringtie_out), file(stringtie_out_gtf)

        output:
        tuple file("${id}.gene.readCounts"), file("${id}.transcripts.readCounts")

        script:
        """
        echo $id $stringtie_out_gtf > prepDE_inp.${id}.txt
        python3 $baseDir/prepDE.py -i prepDE_inp.${id}.txt -g ${id}.gene.readCounts -t ${id}.transcripts.readCounts
        """
}

process stringtie_out_parsing{

        input:
        tuple val(id), file(stringtie_out), file(stringtie_out_gtf)
        tuple file(gene_readCounts), file(transcripts_readCounts)

        output:
        publishDir "${params.outdir}"
        file("${id}.stringtie.out")

        """
        #!/usr/bin/env python2
        import sys
        import subprocess
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pandas"]) 

        import pandas as pd
        
        df_rc_raw = pd.read_csv("${gene_readCounts}",header=0,names=['gene','readCounts'])
        df_rc_raw['gene_ensembl'] = df_rc_raw['gene'].apply(lambda x:x.split('|')[0])
        df_rc_raw['gene_symbol'] = df_rc_raw['gene'].apply(lambda x:x.split('|')[1])
        df_rc = df_rc_raw.loc[:,['gene_ensembl','gene_symbol','readCounts']].set_index('gene_ensembl')

        df_norm_raw = pd.read_csv("${stringtie_out}",sep='\t')
        df_norm = df_norm_raw.loc[:,['Gene ID', 'FPKM', 'TPM']].set_index('Gene ID')

        df_final = df_rc.join(df_norm).reset_index()
        df_final.to_csv("${id}.stringtie.out",sep='\t',index=False,header=True)
        """
}

// this provides an entry point for this main script, so it can be run directly without clone the repo
// using this command: nextflow run <git_acc>/<repo>/<pkg_name>/<main_script>.nf -r <pkg_name>.v<pkg_version> --params-file xxx
workflow {
  stringtie(inp_bam_ch, params.annotation)
  stringtie_rc(stringtie.out)
  stringtie_out_parsing(stringtie.out, stringtie_rc.out)
}
