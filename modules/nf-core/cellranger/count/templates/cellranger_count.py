#!/usr/bin/env python3
"""
Download and rename files for input into cellranger count.
Handles S3, GCS, Azure, HTTP, and local file paths.

Copyright (c) Gregor Sturm 2023 - MIT License
"""

from subprocess import run
from pathlib import Path
from textwrap import dedent
from urllib.parse import urlparse
import shlex
import re
import os


def chunk_iter(seq, size):
    """iterate over `seq` in chunks of `size`"""
    return (seq[pos : pos + size] for pos in range(0, len(seq), size))


def download_file(src: str, dest: Path) -> None:
    """Download a file from various sources (S3, GCS, Azure, HTTP, or local)."""
    parsed = urlparse(src)
    
    if parsed.scheme == 's3':
        run(['aws', 's3', 'cp', src, str(dest), '--only-show-errors'], check=True)
    elif parsed.scheme == 'gs':
        run(['gsutil', 'cp', src, str(dest)], check=True)
    elif parsed.scheme in ('az', 'https') and '.blob.core.windows.net' in src:
        run(['azcopy', 'copy', src, str(dest)], check=True)
    elif parsed.scheme in ('http', 'https'):
        run(['curl', '-sL', src, '-o', str(dest)], check=True)
    elif os.path.isfile(src):
        os.symlink(src, dest)
    else:
        raise ValueError(f"Cannot access file: {src}")


sample_id = "${meta.id}"

# Get reads from the val input (passed as space-separated quoted strings)
reads_str = """${reads_str}"""
# Parse the reads - they come as space-separated quoted paths
import shlex as shlex_parse
read_paths = shlex_parse.split(reads_str)

assert len(read_paths) % 2 == 0, f"Expected even number of reads (R1/R2 pairs), got {len(read_paths)}"

# target directory for downloaded fastqs
fastq_download = Path("./fastq_download")
fastq_download.mkdir(exist_ok=True)

# Download all files first
downloaded_files = []
for i, read_path in enumerate(read_paths):
    filename = Path(read_path).name
    dest = fastq_download / f"{i:03d}_{filename}"
    print(f"Downloading {read_path} to {dest}")
    download_file(read_path, dest)
    downloaded_files.append(dest)

# target directory in which the renamed fastqs will be placed
fastq_all = Path("./fastq_all")
fastq_all.mkdir(exist_ok=True)

# Match R1 in the filename, but only if it is followed by a non-digit or non-character
# match "file_R1.fastq.gz", "file.R1_000.fastq.gz", etc. but
# do not match "SRR12345", "file_INFIXR12", etc
filename_pattern = r"([^a-zA-Z0-9])R1([^a-zA-Z0-9])"

for i, (r1, r2) in enumerate(chunk_iter(downloaded_files, 2), start=1):
    # double escapes are required because nextflow processes this python 'template'
    if re.sub(filename_pattern, r"\\1R2\\2", r1.name.split('_', 1)[1]) != r2.name.split('_', 1)[1]:
        raise AssertionError(
            dedent(
                f"""\
                We expect R1 and R2 of the same sample to have the same filename except for R1/R2.
                This has been checked by replacing "R1" with "R2" in the first filename and comparing it to the second filename.
                If you believe this check shouldn't have failed on your filenames, please report an issue on GitHub!

                Files involved:
                    - {r1}
                    - {r2}
                """
            )
        )
    r1.rename(fastq_all / f"{sample_id}_S1_L{i:03d}_R1_001.fastq.gz")
    r2.rename(fastq_all / f"{sample_id}_S1_L{i:03d}_R2_001.fastq.gz")

# fmt: off
run(
    [
        "cellranger", "count",
        "--id", "${prefix}",
        "--fastqs", str(fastq_all),
        "--transcriptome", "${reference.name}",
        "--localcores", "${task.cpus}",
        "--localmem", "${task.memory.toGiga()}",
        *shlex.split("""${args}"""),
    ],
    check=True,
)
# fmt: on

# Output version information
version = run(
    ["cellranger", "-V"],
    text=True,
    check=True,
    capture_output=True,
).stdout.replace("cellranger cellranger-", "")

# alas, no `pyyaml` pre-installed in the cellranger container
with open("versions.yml", "w") as f:
    f.write('"${task.process}":\\n')
    f.write(f'  cellranger: "{version}"\\n')
