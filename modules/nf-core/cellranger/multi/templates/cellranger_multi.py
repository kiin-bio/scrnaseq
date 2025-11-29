#!/usr/bin/env python3
"""
Download and rename files for input into cellranger multi and run it.
Handles S3, GCS, Azure, HTTP, and local file paths.

Copyright (c) Felipe Almeida 2024 - MIT License
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


sample_id = "${prefix}"
EMPTY_FILE = "EMPTY"

# target directory in which the renamed fastqs will be placed
fastq_all = Path("./fastq_all")
fastq_all.mkdir(exist_ok=True)

# Match R1 in the filename, but only if it is followed by a non-digit or non-character
# match "file_R1.fastq.gz", "file.R1_000.fastq.gz", etc. but
# do not match "SRR12345", "file_INFIXR12", etc
filename_pattern = r"([^a-zA-Z0-9])R1([^a-zA-Z0-9])"

# Fastq paths from val inputs (as space-separated quoted strings)
fastq_inputs = {
    "gex": """${gex_fastqs_str}""",
    "vdj": """${vdj_fastqs_str}""",
    "ab": """${ab_fastqs_str}""",
    "beam": """${beam_fastqs_str}""",
    "cmo": """${cmo_fastqs_str}""",
    "crispr": """${crispr_fastqs_str}""",
}

for modality, fastqs_str in fastq_inputs.items():
    # Parse the fastq paths from the string
    read_paths = shlex.split(fastqs_str) if fastqs_str.strip() else []
    
    if not read_paths:
        continue
        
    # Skip if only contains EMPTY placeholder
    if all(EMPTY_FILE in p for p in read_paths):
        continue
    
    assert len(read_paths) % 2 == 0, f"Expected even number of reads for {modality}, got {len(read_paths)}"

    # Download directory for this modality
    download_dir = Path(f"./fastq_download/{modality}")
    download_dir.mkdir(parents=True, exist_ok=True)
    
    # Download all files
    downloaded_files = []
    for i, read_path in enumerate(read_paths):
        filename = Path(read_path).name
        dest = download_dir / f"{i:03d}_{filename}"
        print(f"Downloading {read_path} to {dest}")
        download_file(read_path, dest)
        downloaded_files.append(dest)

    # target directory in which the renamed fastqs will be placed
    final_dir = fastq_all / modality
    final_dir.mkdir(exist_ok=True)

    for i, (r1, r2) in enumerate(chunk_iter(downloaded_files, 2), start=1):
        # will we rename the files or just move it with the same name to
        # the 'fastq_all' directory which is where files are expected?
        if "${skip_renaming}" == "true":  # nf variables are true/false, which are different from Python
            # Remove the index prefix we added during download
            resolved_name_r1 = r1.name.split('_', 1)[1] if '_' in r1.name else r1.name
            resolved_name_r2 = r2.name.split('_', 1)[1] if '_' in r2.name else r2.name

        else:
            # Get original names without index prefix
            r1_orig = r1.name.split('_', 1)[1] if '_' in r1.name else r1.name
            r2_orig = r2.name.split('_', 1)[1] if '_' in r2.name else r2.name
            
            # double escapes are required because nextflow processes this python 'template'
            if re.sub(filename_pattern, r"\\1R2\\2", r1_orig) != r2_orig:
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

            resolved_name_r1 = f"{sample_id}_S1_L{i:03d}_R1_001.fastq.gz"
            resolved_name_r2 = f"{sample_id}_S1_L{i:03d}_R2_001.fastq.gz"

        # rename or just move
        r1.rename(final_dir / resolved_name_r1)
        r2.rename(final_dir / resolved_name_r2)

#
# fix relative paths from main.nf
#
work_dir = str(Path.cwd()) + "/"
gex_reference_path = "${gex_reference_path}".replace("./", work_dir)
frna_probeset = "${frna_probeset}".replace("./", work_dir)
target_panel = "${target_panel}".replace("./", work_dir)
cmo_reference_path = "${cmo_reference_path}".replace("./", work_dir)
cmo_barcode_path = "${cmo_barcode_path}".replace("./", work_dir)
fb_reference_path = "${fb_reference_path}".replace("./", work_dir)
vdj_reference_path = "${vdj_reference_path}".replace("./", work_dir)
primer_index = "${primer_index}".replace("./", work_dir)
fastq_gex = "${fastq_gex}".replace("./", work_dir)
fastq_vdj = "${fastq_vdj}".replace("./", work_dir)
fastq_antibody = "${fastq_antibody}".replace("./", work_dir)
fastq_beam = "${fastq_beam}".replace("./", work_dir)
fastq_crispr = "${fastq_crispr}".replace("./", work_dir)
fastq_cmo = "${fastq_cmo}".replace("./", work_dir)

#
# generate config file for cellranger multi
#
config_txt = f"""${include_gex}
{gex_reference_path}
{frna_probeset}
${gex_options_filter_probes}
${gex_options_r1_length}
${gex_options_r2_length}
${gex_options_chemistry}
${gex_options_expect_cells}
${gex_options_force_cells}
${gex_options_no_secondary}
${gex_options_no_bam}
${gex_options_check_library_compatibility}
{target_panel}
${gex_options_no_target_umi_filter}
${gex_options_include_introns}
${cmo_options_min_assignment_confidence}
{cmo_reference_path}
{cmo_barcode_path}

${include_fb}
{fb_reference_path}
${fb_options_r1_length}
${fb_options_r2_length}

${include_vdj}
{vdj_reference_path}
{primer_index}
${vdj_options_r1_length}
${vdj_options_r2_length}

[libraries]
fastq_id,fastqs,lanes,feature_types
{fastq_gex}
{fastq_vdj}
{fastq_antibody}
{fastq_beam}
{fastq_crispr}
{fastq_cmo}
"""

#
# check the extra data that is included
#
if len("${include_cmo}") > 0:
    with open("${cmo_csv_text}", "r") as input_conf:
        config_txt = config_txt + "\\n${include_cmo}\\n" + input_conf.read() + "\\n"

if len("${include_beam}") > 0:
    with open("${beam_csv_text}", "r") as input_conf, open("${beam_antigen_csv}", "r") as input_csv:
        config_txt = config_txt + "\\n${include_beam}\\n" + input_conf.read() + "\\n"
        config_txt = config_txt + "[feature]\\n" + input_csv.read() + "\\n"

if len("${include_frna}") > 0:
    with open("${frna_csv_text}", "r") as input_conf:
        config_txt = config_txt + "\\n${include_frna}\\n" + input_conf.read() + "\\n"

# Remove blank lines from config text
config_txt = "\\n".join([line for line in config_txt.split("\\n") if line.strip() != ""])

# Save config file
with open("${config}", "w") as file:
    file.write(config_txt)

#
# run cellranger multi
#
run(
    # fmt: off
    [
        "cellranger",
        "multi",
        f"--id={sample_id}",
        "--csv=${config}",
        "--localcores=${task.cpus}",
        "--localmem=${task.memory.toGiga()}",
        *shlex.split("""${args}"""),
    ],
    # fmt: on
    check=True,
)

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
