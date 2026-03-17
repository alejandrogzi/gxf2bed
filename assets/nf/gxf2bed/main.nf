process GXF2BED {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        '' : 
        'ghcr.io/alejandrogzi/gxf2bed:latest' }"

    input:
    tuple val(meta), path(gxf)

    output:
    tuple val(meta), path("*.bed"),   emit: bed
    path "versions.yml"           ,   emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: gxf.baseName
    """
    gxf2bed \\
        $args \\
        $gxf \\
        ${prefix}.bed
        
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gxf2bed: \$(gxf2bed --version | sed -e "s/gxf2bed v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bed

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gxf2bed: \$(gxf2bed --version | sed -e "s/gxf2bed v//g")
    END_VERSIONS
    """
}
