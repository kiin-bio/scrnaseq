process FASTQC {
    tag "${meta.id}"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/fastqc:0.12.1--hdfd78af_0' :
        'biocontainers/fastqc:0.12.1--hdfd78af_0' }"

    input:
    tuple val(meta), val(reads)

    output:
    tuple val(meta), path("*.html"), emit: html
    tuple val(meta), path("*.zip") , emit: zip
    path  "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args          = task.ext.args ?: ''
    def prefix        = task.ext.prefix ?: "${meta.id}"
    // Convert reads to list if single file
    def reads_list    = reads instanceof List ? reads : [reads]
    def reads_names   = reads_list.withIndex().collect { read, index -> 
        def basename = read.toString().split('/')[-1]
        def extension = basename.endsWith('.gz') ? 'fastq.gz' : 'fastq'
        reads_list.size() == 1 ? "${prefix}.${extension}" : "${prefix}_${index + 1}.${extension}"
    }
    def renamed_files = reads_names.join(' ')
    def reads_str     = reads_list.collect { "\"${it}\"" }.join(' ')

    // The total amount of allocated RAM by FastQC is equal to the number of threads defined (--threads) time the amount of RAM defined (--memory)
    // https://github.com/s-andrews/FastQC/blob/1faeea0412093224d7f6a07f777fad60a5650795/fastqc#L211-L222
    // Dividing the task.memory by task.cpu allows to stick to requested amount of RAM in the label
    def memory_in_mb = task.memory ? task.memory.toUnit('MB').toFloat() / task.cpus : null
    // FastQC memory value allowed range (100 - 10000)
    def fastqc_memory = memory_in_mb > 10000 ? 10000 : (memory_in_mb < 100 ? 100 : memory_in_mb)

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

    # Download reads to local directory with renamed filenames
    reads_arr=(${reads_str})
    names_arr=(${renamed_files})
    for i in "\${!reads_arr[@]}"; do
        download_file "\${reads_arr[\$i]}" "\${names_arr[\$i]}"
    done

    fastqc \\
        ${args} \\
        --threads ${task.cpus} \\
        --memory ${fastqc_memory} \\
        ${renamed_files}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.html
    touch ${prefix}.zip

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        fastqc: \$( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
    END_VERSIONS
    """
}
