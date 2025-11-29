process DOWNLOAD_FASTQS {
    tag "$meta.id"
    label 'process_low'

    // Use amazonlinux with AWS CLI v2 pre-installed - has normal shell entrypoint
    // Unlike amazon/aws-cli which has 'aws' as entrypoint and can't run bash scripts
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/awscli:1.29.10--pyh7cba7a3_1' :
        '686255957239.dkr.ecr.eu-west-1.amazonaws.com/downloadfastq:latest' }"

    env AWS_ACCESS_KEY_ID, "$AWS_ACCESS_KEY_ID" 
    env AWS_SECRET_ACCESS_KEY, "$AWS_SECRET_ACCESS_KEY"
    env AWS_DEFAULT_REGION, 'eu-west-1'

    input:
    tuple val(meta), val(reads)

    output:
    tuple val(meta), path("fastqs/*.fastq.gz"), emit: fastqs
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def reads_list = reads instanceof List ? reads : [reads]
    def reads_str  = reads_list.collect { "\"${it}\"" }.join(' ')
    """
    mkdir -p fastqs

    # Helper function to download files based on URL type
    download_file() {
        local src="\$1"
        local dest="\$2"
        
        if [[ "\$src" == s3://* ]]; then
            aws s3 cp "\$src" "\$dest" --only-show-errors
        elif [[ "\$src" == az://* ]] || [[ "\$src" == https://*.blob.core.windows.net/* ]]; then
            azcopy copy "\$src" "\$dest" 2>/dev/null || az storage blob download --blob-url "\$src" --file "\$dest"
        elif [[ "\$src" == gs://* ]]; then
            gsutil -q cp "\$src" "\$dest"
        elif [[ "\$src" == http://* ]] || [[ "\$src" == https://* ]]; then
            curl -sL "\$src" -o "\$dest"
        elif [[ -f "\$src" ]]; then
            ln -s "\$src" "\$dest"
        else
            echo "ERROR: Cannot access file: \$src" >&2
            exit 1
        fi
    }

    # Download each file
    for read_path in ${reads_str}; do
        filename=\$(basename "\$read_path")
        echo "Downloading \$read_path -> fastqs/\$filename"
        download_file "\$read_path" "fastqs/\$filename"
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aws-cli: \$(aws --version 2>&1 | cut -d/ -f2 | cut -d' ' -f1 || echo "N/A")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p fastqs
    touch fastqs/stub.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        aws-cli: \$(aws --version 2>&1 | cut -d/ -f2 | cut -d' ' -f1 || echo "N/A")
    END_VERSIONS
    """
}

