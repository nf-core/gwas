#!/usr/bin/env bash

set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)
profile=docker
arguments=("$@")

for ((index = 0; index < ${#arguments[@]}; index++)); do
    argument=${arguments[$index]}
    case "$argument" in
        --profile=*) profile=${argument#--profile=} ;;
        --profile)
            if ((index + 1 < ${#arguments[@]})); then
                profile=${arguments[$((index + 1))]}
            fi
            ;;
    esac
done
profile=${profile#+}

fixture_root=$(
    "$script_dir/materialize.sh" --profile "$profile"
)
export GWAS_TEST_FIXTURES=$fixture_root
exec nf-test "$@"
