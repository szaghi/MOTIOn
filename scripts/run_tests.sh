#!/usr/bin/env bash
# Run the MOTIOn test suite.
# Each motion_* binary exits 0 on success, non-zero on failure.
# All tests are MPI-compiled; override the process count with NP=<n>.

NP=${NP:-2}

if [[ -t 1 ]]; then
  RED=$'\033[0;31m'; GREEN=$'\033[0;32m'; BOLD=$'\033[1m'; RESET=$'\033[0m'
else
  RED=''; GREEN=''; BOLD=''; RESET=''
fi

pass=0; fail=0
tmpout=$(mktemp)
trap 'rm -f "$tmpout"' EXIT

for exe in exe/motion_*; do
  [[ -x "$exe" ]] || continue
  name=$(basename "$exe")
  if mpirun -np "$NP" "$exe" > "$tmpout" 2>&1; then
    printf "  ${GREEN}PASS${RESET}  %s\n" "$name"
    pass=$((pass + 1))
  else
    printf "  ${RED}FAIL${RESET}  %s\n" "$name"
    sed 's/^/       /' "$tmpout"
    fail=$((fail + 1))
  fi
done

total=$((pass + fail))
printf "\n${BOLD}%d/%d passed${RESET}\n" "$pass" "$total"
exit $fail
