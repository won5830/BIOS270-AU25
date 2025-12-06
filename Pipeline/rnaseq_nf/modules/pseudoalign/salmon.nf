process SALMON {
    tag "$sample"
    publishDir "${params.outdir}/${sample}/salmon_outs", mode: 'copy', overwrite: true

    input:
      tuple val(sample), path(r1), path(r2), val(condition)
      path index
    output:
      tuple val(sample), path('salmon_outs/quant.sf'), val(condition)
    script:
      """
      set -euo pipefail
      mkdir -p salmon_outs
      salmon quant -i "${index}" -l A \
                   -1 "$r1" -2 "$r2" \
                   -p ${task.cpus} --validateMappings \
                   -o salmon_outs
      """
}


process SALMON_INDEX {

    tag "salmon_index"

    input:
      path transcriptome_fa   

    output:
      path 'salmon_index'     

    script:
      """
      set -euo pipefail

      salmon index \
        -t "${transcriptome_fa}" \
        -i salmon_index \
        -p ${task.cpus}
      """
}

