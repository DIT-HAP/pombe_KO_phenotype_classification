#!/usr/bin/env bash
# =============================================================================
# run_pipeline.sh — Hayles 2013 phenotype classification pipeline
# =============================================================================
#
# Runs the numbered steps in dependency order, not numeric order:
#
#   00 -> 01 -> 03 -> 04 -> 05 -> {02, 06, 07}
#
# 00-05 use the "bioinformatics" conda env; 06-07 need the "cnsplots" env.
# Most steps resolve their relative default paths against the current
# working directory, so this runner always cd's to the repository root.
#
# Every path is passed explicitly, so this file is also the reference for
# which parameters each step uses. Run ./run_pipeline.sh --help for options.
#
# Exit codes: 0 success, 1 a step failed, 2 usage or precondition error.
# =============================================================================

set -euo pipefail

# =============================================================================
# PROJECT ROOT
# =============================================================================

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT_DIR"

# =============================================================================
# CONFIGURATION
# =============================================================================

MAMBA_BIN="${MAMBA_BIN:-mamba}"
BIO_ENV="${BIO_ENV:-bioinformatics}"
PLOT_ENV="${PLOT_ENV:-cnsplots}"

DEFAULT_RELEASE="2026-06-01"
RELEASE="$DEFAULT_RELEASE"

REF_DIR="data/references"
RAW_DIR="data/raw"
RESULTS_DIR="results"

KIM_2010="${REF_DIR}/kim2010/kim2010_supp_table1.xls"
HAYLES_2013="${REF_DIR}/hayles2013/hayles2013_supp_table3.xlsx"
MANUAL_INCONSISTENT="${REF_DIR}/previous_manual_check_of_insistent_phenotypes/Inconsistent_phenotypes_at_25_32_manual.xlsx"
DIT_HAP="${RAW_DIR}/gene_level_fitting_statistics_minusDR.tsv"
GRNA="${RAW_DIR}/260127-all_genes_order1_gRNA_HDdata_fitted_parameters.tsv"
CATEGORY_REVISIONS="data/4_categorized_genes/category_revisions.json"

KIM_PREMERGE="data/1_formatted/Kim2010_paralog_premerge.xlsx"
HAYLES_PREMERGE="data/1_formatted/Hayles2013_phenotypes_premerge.xlsx"
ANNOTATED="data/1_formatted/Hayles2013_phenotypes_annotated.xlsx"
GROUPED="data/3_grouped_genes/Hayles_2013_OB_grouped_genes.xlsx"
CATEGORIZED="data/4_categorized_genes/Hayles_2013_OB_categorized_phenotypes.xlsx"
MERGED="data/5_merged_categories/Hayles_2013_OB_merged_categories.xlsx"
RESULTS_MERGED="${RESULTS_DIR}/Hayles_2013_OB_merged_categories.xlsx"
KEYWORD_PROFILE="data/2_keyword_profile/growth_tier_keyword_profile.xlsx"

# Set from --release in main(); used by steps 00 and 01.
ANNOTATION=""

STEP_ORDER=(00 01 03 04 05 02 06 07)

FROM_STEP=""
TO_STEP=""
ONLY_LIST=""
SKIP_LIST=""
DRY_RUN=0
VERBOSE=0
FORCE=0
LIST_ONLY=0

# =============================================================================
# HELPERS
# =============================================================================

die() {
    printf '[ERROR] %s\n' "$*" >&2
    exit 2
}

usage() {
    cat <<'EOF'
Usage: ./run_pipeline.sh [options]

Run the Hayles 2013 phenotype classification pipeline. The runner cd's to
its own directory (the repository root). Steps run in dependency order, not
numeric order:

    00 -> 01 -> 03 -> 04 -> 05 -> {02, 06, 07}

00-05 use BIO_ENV (default: bioinformatics); 06-07 use PLOT_ENV
(default: cnsplots).

Selection (applied in execution order, not numeric order):
  --from STEP       Start at STEP (e.g. --from 04)
  --to STEP         Stop after STEP (e.g. --to 05)
  --only LIST       Run only these comma-separated steps (e.g. --only 05,06)
  --skip LIST       Skip these comma-separated steps (e.g. --skip 06,07)
  --list            Print the step table and exit

Execution:
  --dry-run         Print every command and path without running it
  --verbose         Pass --verbose to every step
  --force           Pass --force to step 00 (re-download PomBase)
  --release DATE    PomBase release (default: 2026-06-01); also sets the
                    annotation path passed to step 01
  -h, --help        Show this help

Environment variables:
  MAMBA_BIN   mamba executable (default: mamba)
  BIO_ENV     conda env for 00-05 (default: bioinformatics)
  PLOT_ENV    conda env for 06-07 (default: cnsplots)

Examples:
  ./run_pipeline.sh
  ./run_pipeline.sh --from 04
  ./run_pipeline.sh --to 05
  ./run_pipeline.sh --only 05,06
  ./run_pipeline.sh --skip 06,07
  ./run_pipeline.sh --dry-run --only 01
  ./run_pipeline.sh --release 2026-05-01 --force --only 00,01

Selections do NOT run prerequisites automatically. Each selected step checks
its inputs before it starts; the error names the step that produces the file.
EOF
}

normalize_step() {
    case "$1" in
        0|00) printf '00' ;;
        1|01) printf '01' ;;
        2|02) printf '02' ;;
        3|03) printf '03' ;;
        4|04) printf '04' ;;
        5|05) printf '05' ;;
        6|06) printf '06' ;;
        7|07) printf '07' ;;
        *) return 1 ;;
    esac
}

step_index() {
    local step=$1 i
    for ((i = 0; i < ${#STEP_ORDER[@]}; i++)); do
        if [[ "${STEP_ORDER[$i]}" == "$step" ]]; then
            printf '%s' "$i"
            return 0
        fi
    done
    return 1
}

step_env() {
    case "$1" in
        06|07) printf '%s' "$PLOT_ENV" ;;
        *) printf '%s' "$BIO_ENV" ;;
    esac
}

step_script() {
    case "$1" in
        00) printf '%s' "scripts/00_download_pombase_annotation.py" ;;
        01) printf '%s' "scripts/01_format_hayles_phenotypes.py" ;;
        02) printf '%s' "scripts/02_keyword_profile.py" ;;
        03) printf '%s' "scripts/03_group_genes.py" ;;
        04) printf '%s' "scripts/04_categorize_phenotypes.py" ;;
        05) printf '%s' "scripts/05_merge_categories.py" ;;
        06) printf '%s' "scripts/06_plot_dr_distribution.py" ;;
        07) printf '%s' "scripts/07_plot_dr_distribution_pure.py" ;;
    esac
}

step_description() {
    case "$1" in
        00) printf '%s' "Download the PomBase annotation" ;;
        01) printf '%s' "Format Hayles phenotypes (Kim paralogs)" ;;
        02) printf '%s' "Profile keywords per Growth_tier" ;;
        03) printf '%s' "Group genes by temperature consistency" ;;
        04) printf '%s' "Categorise growth phenotypes" ;;
        05) printf '%s' "Merge categories and re-rank tiers" ;;
        06) printf '%s' "Plot DR/um distributions (original + revised)" ;;
        07) printf '%s' "Plot the five pure categories" ;;
    esac
}

format_duration() {
    local seconds=${1:-0}
    if ((seconds < 60)); then
        printf '%ds' "$seconds"
    elif ((seconds < 3600)); then
        printf '%dm %02ds' "$((seconds / 60))" "$((seconds % 60))"
    else
        printf '%dh %02dm %02ds' \
            "$((seconds / 3600))" \
            "$(((seconds % 3600) / 60))" \
            "$((seconds % 60))"
    fi
}

print_command() {
    local arg first=1
    for arg in "$@"; do
        if ((first)); then
            first=0
        else
            printf ' '
        fi
        printf '%q' "$arg"
    done
    printf '\n'
}

build_step_set() {
    # $1 = comma-separated list; sets STEP_SET to ",00,05,"
    local list=$1 item normalized
    STEP_SET=","
    local IFS=','
    for item in $list; do
        item=${item// /}
        [[ -z "$item" ]] && continue
        if ! normalized=$(normalize_step "$item"); then
            die "unknown step: '$item' (valid: 00 01 02 03 04 05 06 07)"
        fi
        STEP_SET="${STEP_SET}${normalized},"
    done
}

# =============================================================================
# ARGUMENT PARSING
# =============================================================================

parse_args() {
    while [[ $# -gt 0 ]]; do
        case "$1" in
            --from)
                [[ $# -ge 2 && -n "$2" ]] || die "--from needs a step value"
                FROM_STEP=$2
                shift 2
                ;;
            --to)
                [[ $# -ge 2 && -n "$2" ]] || die "--to needs a step value"
                TO_STEP=$2
                shift 2
                ;;
            --only)
                [[ $# -ge 2 && -n "$2" ]] || die "--only needs a non-empty comma-separated list"
                ONLY_LIST=$2
                shift 2
                ;;
            --skip)
                [[ $# -ge 2 && -n "$2" ]] || die "--skip needs a non-empty comma-separated list"
                SKIP_LIST=$2
                shift 2
                ;;
            --release)
                [[ $# -ge 2 ]] || die "--release needs a version (YYYY-MM-DD)"
                RELEASE=$2
                shift 2
                ;;
            --dry-run)
                DRY_RUN=1
                shift
                ;;
            --verbose)
                VERBOSE=1
                shift
                ;;
            --force)
                FORCE=1
                shift
                ;;
            --list)
                LIST_ONLY=1
                shift
                ;;
            -h|--help)
                usage
                exit 0
                ;;
            *)
                die "unknown option: '$1' (try --help)"
                ;;
        esac
    done
}

validate_release() {
    if [[ ! "$RELEASE" =~ ^[0-9]{4}-[0-9]{2}-[0-9]{2}$ ]]; then
        die "invalid --release '$RELEASE' (expected YYYY-MM-DD)"
    fi
}

# =============================================================================
# STEP SELECTION
# =============================================================================

select_steps() {
    local step from_idx to_idx i

    if [[ -n "$ONLY_LIST" ]]; then
        if [[ -n "$FROM_STEP" || -n "$TO_STEP" || -n "$SKIP_LIST" ]]; then
            die "--only cannot be combined with --from, --to, or --skip"
        fi
    fi

    if [[ -n "$FROM_STEP" ]]; then
        local raw=$FROM_STEP
        if ! FROM_STEP=$(normalize_step "$raw"); then
            die "unknown step for --from: '$raw'"
        fi
    fi
    if [[ -n "$TO_STEP" ]]; then
        local raw=$TO_STEP
        if ! TO_STEP=$(normalize_step "$raw"); then
            die "unknown step for --to: '$raw'"
        fi
    fi

    SELECTED=()
    if [[ -n "$ONLY_LIST" ]]; then
        build_step_set "$ONLY_LIST"
        local only_set=$STEP_SET
        for step in "${STEP_ORDER[@]}"; do
            if [[ "$only_set" == *",$step,"* ]]; then
                SELECTED+=("$step")
            fi
        done
    else
        from_idx=0
        to_idx=$((${#STEP_ORDER[@]} - 1))
        if [[ -n "$FROM_STEP" ]]; then
            from_idx=$(step_index "$FROM_STEP")
        fi
        if [[ -n "$TO_STEP" ]]; then
            to_idx=$(step_index "$TO_STEP")
        fi
        if ((from_idx > to_idx)); then
            die "--from $FROM_STEP comes after --to $TO_STEP in execution order"
        fi

        SKIP_SET=""
        if [[ -n "$SKIP_LIST" ]]; then
            build_step_set "$SKIP_LIST"
            SKIP_SET=$STEP_SET
        fi

        for ((i = from_idx; i <= to_idx; i++)); do
            step=${STEP_ORDER[$i]}
            if [[ -n "$SKIP_SET" && "$SKIP_SET" == *",$step,"* ]]; then
                continue
            fi
            SELECTED+=("$step")
        done
    fi

    if ((${#SELECTED[@]} == 0)); then
        die "no steps selected"
    fi
}

# =============================================================================
# PREFLIGHT
# =============================================================================

check_mamba() {
    if ! command -v "$MAMBA_BIN" >/dev/null 2>&1; then
        die "executable not found: $MAMBA_BIN (override with MAMBA_BIN=...)"
    fi
}

check_scripts() {
    local step script
    for step in "${SELECTED[@]}"; do
        script=$(step_script "$step")
        if [[ ! -f "$script" ]]; then
            die "script not found: $script"
        fi
    done
}

require_file() {
    local step=$1 path=$2 hint=$3
    if [[ ! -f "$path" ]]; then
        printf '[ERROR] step %s is missing an input:\n' "$step" >&2
        printf '          %s\n' "$path" >&2
        printf '          %s\n' "$hint" >&2
        return 1
    fi
}

check_step_inputs() {
    local step=$1 missing=0

    case "$step" in
        00)
            ;;
        01)
            require_file "$step" "$KIM_2010" \
                "Kim 2010 Table S1 (external; see AGENTS.md)" || missing=1
            require_file "$step" "$HAYLES_2013" \
                "Hayles 2013 main table (external; see AGENTS.md)" || missing=1
            require_file "$step" "$ANNOTATION" \
                "PomBase annotation; run step 00 or download it manually" || missing=1
            ;;
        02)
            require_file "$step" "$RESULTS_MERGED" \
                "produced by step 05" || missing=1
            ;;
        03)
            require_file "$step" "$ANNOTATED" \
                "produced by step 01" || missing=1
            ;;
        04)
            require_file "$step" "$GROUPED" \
                "produced by step 03" || missing=1
            require_file "$step" "$MANUAL_INCONSISTENT" \
                "hand-curated manual annotations (tracked in git)" || missing=1
            # Step 04 reads DIT_HAP through its DEFAULT_DIT_HAP constant;
            # there is no --dit-hap flag.
            require_file "$step" "$DIT_HAP" \
                "fitted DIT-HAP parameters (data/raw/)" || missing=1
            ;;
        05)
            require_file "$step" "$CATEGORIZED" \
                "produced by step 04" || missing=1
            require_file "$step" "$CATEGORY_REVISIONS" \
                "hand-curated category orders and merges (tracked in git)" || missing=1
            require_file "$step" "$DIT_HAP" \
                "fitted DIT-HAP parameters (data/raw/)" || missing=1
            ;;
        06)
            require_file "$step" "$MERGED" \
                "produced by step 05" || missing=1
            require_file "$step" "$DIT_HAP" \
                "fitted DIT-HAP parameters (data/raw/)" || missing=1
            require_file "$step" "$GRNA" \
                "fitted gRNA parameters (data/raw/)" || missing=1
            require_file "$step" "$CATEGORY_REVISIONS" \
                "hand-curated category orders and merges (tracked in git)" || missing=1
            ;;
        07)
            require_file "$step" "$MERGED" \
                "produced by step 05" || missing=1
            require_file "$step" "$DIT_HAP" \
                "fitted DIT-HAP parameters (data/raw/)" || missing=1
            require_file "$step" "$GRNA" \
                "fitted gRNA parameters (data/raw/)" || missing=1
            require_file "$step" "$CATEGORY_REVISIONS" \
                "hand-curated category orders and merges (tracked in git)" || missing=1
            ;;
    esac

    if ((missing)); then
        die "step $step cannot start; see the missing inputs above"
    fi
}

# =============================================================================
# RUN
# =============================================================================

print_plan() {
    local step

    printf 'Repository: %s\n' "$ROOT_DIR"
    printf 'PomBase:    %s\n' "$RELEASE"
    if ((DRY_RUN)); then
        printf 'Mode:       dry run (nothing executed)\n'
    else
        printf 'Mode:       execute\n'
    fi
    printf 'Selected:   '
    printf '%s ' "${SELECTED[@]}"
    printf '\n\n'

    printf '  %-4s %-16s %s\n' "STEP" "ENV" "DESCRIPTION"
    for step in "${SELECTED[@]}"; do
        printf '  %-4s %-16s %s\n' \
            "$step" "$(step_env "$step")" "$(step_description "$step")"
    done
    printf '\n'

    if ((FORCE)) && [[ " ${SELECTED[*]} " != *" 00 "* ]]; then
        printf '[WARN] --force has no effect because step 00 is not selected.\n\n'
    fi
    if ((DRY_RUN)); then
        printf '[DRY-RUN] Input checks are skipped; commands are printed only.\n\n'
    fi
}

run_step() {
    local step=$1 pos=$2 total=$3
    local env script desc label
    local -a args=()
    local -a cmd=()
    local start rc resume

    env=$(step_env "$step")
    script=$(step_script "$step")
    desc=$(step_description "$step")

    case "$step" in
        00)
            args=(--release "$RELEASE" --outdir "$REF_DIR")
            if ((FORCE)); then
                args+=(--force)
            fi
            ;;
        01)
            args=(--kim2010 "$KIM_2010"
                  --hayles2013 "$HAYLES_2013"
                  --annotation "$ANNOTATION"
                  --kim-premerge-output "$KIM_PREMERGE"
                  --hayles-premerge-output "$HAYLES_PREMERGE"
                  --output "$ANNOTATED")
            ;;
        02)
            args=(--input "$RESULTS_MERGED" --output "$KEYWORD_PROFILE")
            ;;
        03)
            args=(--input "$ANNOTATED" --output "$GROUPED")
            ;;
        04)
            args=(--grouped "$GROUPED"
                  --manual "$MANUAL_INCONSISTENT"
                  --output "$CATEGORIZED")
            ;;
        05)
            args=(--input "$CATEGORIZED"
                  --output "$MERGED"
                  --results "$RESULTS_MERGED"
                  --revisions "$CATEGORY_REVISIONS"
                  --dit-hap "$DIT_HAP")
            ;;
        06)
            args=(--merged "$MERGED"
                  --dit-hap "$DIT_HAP"
                  --grna "$GRNA"
                  --revisions "$CATEGORY_REVISIONS"
                  --output-dir "$RESULTS_DIR")
            ;;
        07)
            args=(--merged "$MERGED"
                  --dit-hap "$DIT_HAP"
                  --grna "$GRNA"
                  --revisions "$CATEGORY_REVISIONS"
                  --output-dir "$RESULTS_DIR")
            ;;
        *)
            die "internal error: unknown step '$step'"
            ;;
    esac

    if ((VERBOSE)); then
        args+=(--verbose)
    fi

    if ((DRY_RUN)); then
        label="DRY-RUN"
    else
        label="STEP"
    fi

    printf '\n[%s %d/%d] %s  %s\n' "$label" "$pos" "$total" "$step" "$desc"
    printf '  env: %s\n' "$env"
    printf '  cmd: '
    cmd=("$MAMBA_BIN" run -n "$env" python "$script" "${args[@]}")
    print_command "${cmd[@]}"

    if ((DRY_RUN)); then
        return 0
    fi

    start=$SECONDS
    rc=0
    "${cmd[@]}" || rc=$?

    if ((rc != 0)); then
        if ((rc == 130)); then
            printf '\n[ABORT] step %s interrupted\n' "$step" >&2
            exit 130
        fi
        resume="./run_pipeline.sh --from $step"
        if [[ "$RELEASE" != "$DEFAULT_RELEASE" ]]; then
            resume="$resume --release $RELEASE"
        fi
        printf '\n[FAIL] step %s exited with status %d\n' "$step" "$rc" >&2
        printf '       resume with: %s\n' "$resume" >&2
        return 1
    fi

    printf '[DONE] step %s in %s\n' \
        "$step" "$(format_duration $((SECONDS - start)))"
}

run_selected() {
    local total=${#SELECTED[@]}
    local pos=0 step
    local pipeline_start=$SECONDS

    for step in "${SELECTED[@]}"; do
        pos=$((pos + 1))
        if ((DRY_RUN == 0)); then
            check_step_inputs "$step"
        fi
        if ! run_step "$step" "$pos" "$total"; then
            exit 1
        fi
    done

    if ((DRY_RUN)); then
        printf '\n[OK] %d command(s) printed (dry run)\n' "$total"
    else
        printf '\n[OK] %d step(s) completed in %s\n' \
            "$total" "$(format_duration $((SECONDS - pipeline_start)))"
    fi
}

list_steps() {
    local step
    printf 'Execution order (02 runs after 05, despite its number):\n\n'
    printf '  %-4s %-16s %s\n' "STEP" "ENV" "DESCRIPTION"
    for step in "${STEP_ORDER[@]}"; do
        printf '  %-4s %-16s %s\n' \
            "$step" "$(step_env "$step")" "$(step_description "$step")"
        printf '       %s\n' "$(step_script "$step")"
    done
    printf '\nUse --dry-run to print every command and path.\n'
}

# =============================================================================
# MAIN
# =============================================================================

main() {
    parse_args "$@"
    validate_release
    ANNOTATION="${REF_DIR}/pombase-${RELEASE}_gene_IDs_names_products.tsv"

    if ((LIST_ONLY)); then
        list_steps
        return 0
    fi

    select_steps
    check_scripts
    if ((DRY_RUN == 0)); then
        check_mamba
    fi
    print_plan
    run_selected
}

main "$@"
