#!/usr/bin/env bash

for grm_extension in .grm.bin .grm.N.bin .grm.id; do
    cat \$(ls *"\${grm_extension}" | sort -V) > "${prefix}\${grm_extension}.tmp"
    mv "${prefix}\${grm_extension}.tmp" "${prefix}\${grm_extension}"
done

cat <<-END_VERSIONS > versions.yml
"${task.process}":
    coreutils: \$(cat --version | head -n 1 | cut -d ' ' -f 4)
END_VERSIONS
