#!/bin/bash
# MOCCa pre-commit hook
# This script runs fortitude (linter) and fprettify (formatter check)
# on staged Fortran files before allowing a commit.
#
# Usage:
#   1. Make this script executable: chmod +x scripts/pre-commit-hook.sh
#   2. Create a symlink in .git/hooks/:
#      ln -sf ../../scripts/pre-commit-hook.sh .git/hooks/pre-commit
#
# Or to set it up automatically:
#   ./scripts/pre-commit-hook.sh --install

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Check if this is an install request
if [[ "${1:-}" == "--install" ]]; then
    echo "Installing pre-commit hook..."
    mkdir -p .git/hooks
    ln -sf ../../scripts/pre-commit-hook.sh .git/hooks/pre-commit
    chmod +x .git/hooks/pre-commit
    echo "✓ Pre-commit hook installed successfully!"
    echo "  Run 'git commit' to test it."
    exit 0
fi

# Check if this is an uninstall request
if [[ "${1:-}" == "--uninstall" ]]; then
    echo "Uninstalling pre-commit hook..."
    rm -f .git/hooks/pre-commit
    echo "✓ Pre-commit hook uninstalled."
    exit 0
fi

# Main hook logic
echo -e "${YELLOW}Running MOCCa pre-commit checks...${NC}"

# Get list of staged Fortran files
STAGED_FILES=$(git diff --cached --name-only --diff-filter=ACM | grep -E '\.(f90|f|f95|for)$' || true)

if [[ -z "$STAGED_FILES" ]]; then
    echo -e "${GREEN}No Fortran files staged for commit.${NC}"
    exit 0
fi

echo -e "${YELLOW}Checking ${#STAGED_FILES} Fortran file(s)...${NC}"

# Temporary file for collecting errors
ERRORS_FILE=$(mktemp)
trap "rm -f $ERRORS_FILE" EXIT

# Run fortitude on staged files
echo -e "\n${YELLOW}=== Running fortitude linter ===${NC}"
if fortitude --config-file .fortitude.toml check $STAGED_FILES 2>&1 | tee -a "$ERRORS_FILE"; then
    echo -e "${GREEN}✓ Fortitude passed${NC}"
else
    echo -e "${RED}✗ Fortitude found issues (see above)${NC}" >&2
    HAS_ERRORS=1
fi

# Run fprettify check on staged files
echo -e "\n${YELLOW}=== Running fprettify format check ===${NC}"
FPRETTIFY_OUTPUT=$(mktemp)
fprettify -d -c .fprettify.rc $STAGED_FILES > "$FPRETTIFY_OUTPUT" 2>&1 || true

if [[ ! -s "$FPRETTIFY_OUTPUT" ]]; then
    echo -e "${GREEN}✓ All files are properly formatted${NC}"
else
    echo -e "${RED}✗ Files need formatting${NC}" >&2
    echo "Formatting diff:"
    cat "$FPRETTIFY_OUTPUT"
    echo -e "\n${YELLOW}To fix formatting, run:${NC}"
    echo "  fprettify -r -c .fprettify.rc $(echo $STAGED_FILES | tr ' ' '\n' | sed 's/^/  /')"
    HAS_ERRORS=1
fi

rm -f "$FPRETTIFY_OUTPUT"

# Check if there were any errors
if [[ -s "$ERRORS_FILE" && ${HAS_ERRORS:-0} -eq 1 ]]; then
    echo -e "\n${RED}Pre-commit checks failed!${NC}" >&2
    echo -e "Please fix the issues above and try committing again.\n"
    exit 1
fi

echo -e "\n${GREEN}All pre-commit checks passed!${NC}"
exit 0
