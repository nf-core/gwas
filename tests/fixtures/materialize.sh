#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
project_dir=$(cd "$script_dir/../.." && pwd -P)
profile=docker

usage() {
    printf 'Usage: %s [--profile PROFILE] [--cache-dir DIRECTORY] [--source-root PATH_OR_URL]\n' "$0" >&2
}

cache_dir=
source_root=${GWAS_FIXTURE_SOURCE:-}
while [[ $# -gt 0 ]]; do
    case "$1" in
        --profile)
            [[ $# -ge 2 ]] || { usage; exit 2; }
            profile=$2
            shift 2
            ;;
        --cache-dir)
            [[ $# -ge 2 ]] || { usage; exit 2; }
            cache_dir=$2
            shift 2
            ;;
        --source-root)
            [[ $# -ge 2 ]] || { usage; exit 2; }
            source_root=$2
            shift 2
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            usage
            exit 2
            ;;
    esac
done

worker_root=${NFT_WORKDIR:-"$project_dir/.nf-test"}
case "$worker_root" in
    '~') worker_root=$HOME ;;
    '~/'*) worker_root="$HOME/${worker_root#\~/}" ;;
esac
cache_dir=${cache_dir:-"$worker_root/gwas-fixtures"}
mkdir -p "$cache_dir"
cache_dir=$(cd "$cache_dir" && pwd -P)

canonical_root=${source_root:-https://raw.githubusercontent.com/nf-core/test-datasets/gwas/}
canonical_root=${canonical_root%/}

download_root=$(mktemp -d "${TMPDIR:-/tmp}/gwas-fixture-source.XXXXXX")
run_root=
published_root=
invalid_root=

cleanup() {
    rm -rf -- "$download_root"
    [[ -z "$run_root" ]] || rm -rf -- "$run_root"
    [[ -z "$published_root" ]] || rm -rf -- "$published_root"
    [[ -z "$invalid_root" ]] || rm -rf -- "$invalid_root"
}
trap cleanup EXIT

canonical_files=(
    results/fixtures/genotypes/example_all.vcf.gz
    results/fixtures/pheno_cov/example.pheno
    results/fixtures/pheno_cov/example.qcovar
    results/fixtures/pheno_cov/example.catcovar
)

copy_canonical() {
    local relative_path=$1
    local target=$2
    mkdir -p "$(dirname "$target")"
    if [[ "$canonical_root" == http://* || "$canonical_root" == https://* ]]; then
        curl --fail --location --silent --show-error \
            "$canonical_root/$relative_path" \
            --output "$target"
    else
        local source_path="$canonical_root/$relative_path"
        [[ -f "$source_path" ]] || {
            printf 'Missing canonical GWAS fixture: %s\n' "$source_path" >&2
            exit 1
        }
        cp "$source_path" "$target"
    fi
}

for relative_path in "${canonical_files[@]}"; do
    copy_canonical "$relative_path" "$download_root/$relative_path"
done

sha256_value() {
    if command -v sha256sum >/dev/null 2>&1; then
        sha256sum "$1" | cut -d ' ' -f 1
    else
        shasum -a 256 "$1" | cut -d ' ' -f 1
    fi
}

contract_files=(
    "$script_dir/materialize.sh"
    "$script_dir/materialize.nf"
    "$script_dir/materialize.config"
    "$project_dir/modules/local/plink2/vcf/main.nf"
    "$project_dir/modules/local/plink2/vcf/environment.yml"
    "$project_dir/modules/local/plink2/makebed/main.nf"
    "$project_dir/modules/local/plink2/makebed/environment.yml"
)

digest_input=$(mktemp "${TMPDIR:-/tmp}/gwas-fixture-digest.XXXXXX")
for relative_path in "${canonical_files[@]}"; do
    sha256_value "$download_root/$relative_path" >> "$digest_input"
done
for contract_file in "${contract_files[@]}"; do
    sha256_value "$contract_file" >> "$digest_input"
done
printf '%s\n' "$profile" >> "$digest_input"
digest=$(sha256_value "$digest_input")
rm -f -- "$digest_input"

final_root="$cache_dir/$digest"
manifest_name=.complete.sha256
required_derivatives=(
    results/fixtures/genotypes/example_all.pgen
    results/fixtures/genotypes/example_all.psam
    results/fixtures/genotypes/example_all.pvar
    results/fixtures/genotypes/example_chr1.pgen
    results/fixtures/genotypes/example_chr1.psam
    results/fixtures/genotypes/example_chr1.pvar
    results/fixtures/genotypes/example_all.bed
    results/fixtures/genotypes/example_all.bim
    results/fixtures/genotypes/example_all.fam
)

fixture_files=("${canonical_files[@]}" "${required_derivatives[@]}")

verify_fixture_root() {
    local root=$1
    local manifest="$root/$manifest_name"
    [[ -f "$manifest" ]] || return 1
    [[ $(wc -l < "$manifest") -eq ${#fixture_files[@]} ]] || return 1
    local relative_path
    local expected
    for relative_path in "${fixture_files[@]}"; do
        [[ -f "$root/$relative_path" ]] || return 1
        expected=$(awk -v path="$relative_path" '$2 == path { print $1 }' "$manifest")
        [[ "$expected" =~ ^[[:xdigit:]]{64}$ ]] || return 1
        [[ $(sha256_value "$root/$relative_path") == "$expected" ]] || return 1
    done
}

cache_is_complete() {
    verify_fixture_root "$final_root"
}

if cache_is_complete; then
    printf '%s\n' "$final_root"
    exit 0
fi

command -v flock >/dev/null 2>&1 || {
    printf 'Fixture materialization requires flock for process-lifetime cache locking\n' >&2
    exit 1
}
lock_file="$cache_dir/$digest.lock"
exec 9>"$lock_file"
if ! flock -w 600 9; then
    printf 'Timed out waiting for fixture materialization lock: %s\n' "$lock_file" >&2
    exit 1
fi

# Another writer may have published between the optimistic check and lock acquisition.
if cache_is_complete; then
    printf '%s\n' "$final_root"
    exit 0
fi

published_root=$(mktemp -d "$cache_dir/.$digest.publish.XXXXXX")
run_root=$(mktemp -d "$cache_dir/.$digest.run.XXXXXX")

for relative_path in "${canonical_files[@]}"; do
    mkdir -p "$(dirname "$published_root/$relative_path")"
    cp "$download_root/$relative_path" "$published_root/$relative_path"
done

(
    cd "$project_dir"
    nextflow -log "$run_root/nextflow.log" run "$script_dir/materialize.nf" \
        -c "$script_dir/materialize.config" \
        -profile "$profile" \
        -work-dir "$run_root/work" \
        --fixture_source_vcf "$published_root/results/fixtures/genotypes/example_all.vcf.gz" \
        --fixture_output_dir "$published_root/results/fixtures/genotypes" \
        -ansi-log false
) >&2

for relative_path in "${required_derivatives[@]}"; do
    [[ -f "$published_root/$relative_path" ]] || {
        printf 'Materializer did not publish %s\n' "$relative_path" >&2
        exit 1
    }
done

(
    cd "$published_root"
    for relative_path in "${fixture_files[@]}"; do
        printf '%s  %s\n' "$(sha256_value "$relative_path")" "$relative_path"
    done
) > "$published_root/$manifest_name"
verify_fixture_root "$published_root" || {
    printf 'Materialized fixture manifest verification failed\n' >&2
    exit 1
}

if [[ -e "$final_root" ]]; then
    invalid_root="$cache_dir/.$digest.invalid.$BASHPID"
    mv "$final_root" "$invalid_root"
fi
mv "$published_root" "$final_root"
published_root=
if [[ -n "$invalid_root" ]]; then
    rm -rf -- "$invalid_root"
    invalid_root=
fi

printf '%s\n' "$final_root"
