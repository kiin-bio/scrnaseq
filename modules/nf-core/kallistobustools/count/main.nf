process KALLISTOBUSTOOLS_COUNT {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/kb-python:0.28.2--pyhdfd78af_2' :
        'biocontainers/kb-python:0.28.2--pyhdfd78af_2' }"

    input:
    tuple val(meta), val(reads)
    path  index
    path  t2g
    path  t1c
    path  t2c
    val   technology
    val   workflow_mode

    output:
    tuple val(meta), path ("*.count")   , emit: count
    path "versions.yml"                 , emit: versions
    path "*.count/*/*.mtx"              , emit: matrix //Ensure that kallisto finished and produced outputs

    when:
    task.ext.when == null || task.ext.when

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def cdna       = t1c ? "-c1 $t1c" : ''
    def intron     = t2c ? "-c2 $t2c" : ''
    def memory     = task.memory.toGiga() - 1
    def reads_list = reads instanceof List ? reads : [reads]
    def reads_str  = reads_list.collect { "\"${it}\"" }.join(' ')
    """
    # Helper function to download files based on URL type
    download_file() {
        local src="\$1"
        local dest="\$2"
        
        if [[ "\$src" == s3://* ]]; then
            aws s3 cp "\$src" "\$dest" --only-show-errors
        elif [[ "\$src" == az://* ]] || [[ "\$src" == https://*.blob.core.windows.net/* ]]; then
            azcopy copy "\$src" "\$dest"
        elif [[ "\$src" == gs://* ]]; then
            gsutil cp "\$src" "\$dest"
        elif [[ "\$src" == http://* ]] || [[ "\$src" == https://* ]]; then
            curl -sL "\$src" -o "\$dest"
        elif [[ -f "\$src" ]]; then
            ln -s "\$src" "\$dest"
        else
            echo "ERROR: Cannot access file: \$src" >&2
            exit 1
        fi
    }

    # Download reads to local directory
    mkdir -p fastq_dir
    local_reads=""
    for read_path in ${reads_str}; do
        filename=\$(basename "\$read_path")
        download_file "\$read_path" "fastq_dir/\$filename"
        local_reads="\$local_reads fastq_dir/\$filename"
    done

    kb \\
        count \\
        -t $task.cpus \\
        -i $index \\
        -g $t2g \\
        $cdna \\
        $intron \\
        -x $technology \\
        --workflow $workflow_mode \\
        $args \\
        -o ${prefix}.count \\
        -m ${memory}G \\
        \$local_reads

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    """
    mkdir -p ${prefix}.count/counts_unfiltered/
    touch ${prefix}.count/counts_unfiltered/cells_x_genes.mtx

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        kallistobustools: \$(echo \$(kb --version 2>&1) | sed 's/^.*kb_python //;s/positional arguments.*\$//')
    END_VERSIONS
    """
}
