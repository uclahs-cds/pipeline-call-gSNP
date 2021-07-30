process run_validate {
    container params.docker_image_validate

    input:
    path file_to_validate

    """
    set -euo pipefail
    python -m validate -t file-input ${file_to_validate}
    """
}
