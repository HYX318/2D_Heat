#!/bin/bash

# Script to validate GitHub Actions CI/CD configuration

set -e

echo "=== Validating GitHub Actions CI/CD Configuration ==="
echo ""

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORKFLOWS_DIR="${PROJECT_ROOT}/.github/workflows"

echo "Project root: ${PROJECT_ROOT}"
echo "Workflows directory: ${WORKFLOWS_DIR}"
echo ""

# Check if workflows directory exists
if [ ! -d "$WORKFLOWS_DIR" ]; then
    echo -e "${RED}✗ Workflows directory not found${NC}"
    exit 1
fi
echo -e "${GREEN}✓ Workflows directory exists${NC}"

# Check for required workflow files
REQUIRED_WORKFLOWS=("ci.yml" "benchmark.yml" "codeql.yml")
for workflow in "${REQUIRED_WORKFLOWS[@]}"; do
    if [ -f "${WORKFLOWS_DIR}/${workflow}" ]; then
        echo -e "${GREEN}✓ Found ${workflow}${NC}"
    else
        echo -e "${RED}✗ Missing ${workflow}${NC}"
        exit 1
    fi
done

echo ""
echo "=== Validating YAML Syntax ==="

# Check if yq is available for YAML validation
if command -v yq &> /dev/null; then
    for workflow in {ci,benchmark,codeql}.yml; do
        if yq eval '.' "${WORKFLOWS_DIR}/${workflow}" > /dev/null 2>&1; then
            echo -e "${GREEN}✓ ${workflow} - Valid YAML${NC}"
        else
            echo -e "${RED}✗ ${workflow} - Invalid YAML${NC}"
        fi
    done
else
    echo -e "${YELLOW}⚠ yq not found. Install yq for YAML validation: brew install yq${NC}"
fi

echo ""
echo "=== Checking Workflow Configuration ==="

# Check CI workflow
echo "Checking ci.yml..."
if grep -q "timeout-minutes: 30" "${WORKFLOWS_DIR}/ci.yml"; then
    echo -e "${GREEN}✓ Timeout configured (30 minutes)${NC}"
else
    echo -e "${YELLOW}⚠ Timeout not configured${NC}"
fi

if grep -q "TMPDIR=/tmp" "${WORKFLOWS_DIR}/ci.yml"; then
    echo -e "${GREEN}✓ macOS TMPDIR configured${NC}"
else
    echo -e "${YELLOW}⚠ macOS TMPDIR not configured${NC}"
fi

if grep -q "actions/cache@v3" "${WORKFLOWS_DIR}/ci.yml"; then
    echo -e "${GREEN}✓ Caching enabled${NC}"
else
    echo -e "${YELLOW}⚠ Caching not enabled${NC}"
fi

# Check benchmark workflow
echo "Checking benchmark.yml..."
if grep -q "cron" "${WORKFLOWS_DIR}/benchmark.yml"; then
    echo -e "${GREEN}✓ Scheduled benchmark configured${NC}"
else
    echo -e "${YELLOW}⚠ Scheduled benchmark not configured${NC}"
fi

# Check CodeQL workflow
echo "Checking codeql.yml..."
if grep -q "codeql-action" "${WORKFLOWS_DIR}/codeql.yml"; then
    echo -e "${GREEN}✓ CodeQL action configured${NC}"
else
    echo -e "${YELLOW}⚠ CodeQL action not configured${NC}"
fi

echo ""
echo "=== Checking Matrix Configuration ==="

# Extract OS from CI workflow
if grep -q "os: \[ubuntu-latest, macos-latest\]" "${WORKFLOWS_DIR}/ci.yml"; then
    echo -e "${GREEN}✓ OS matrix: ubuntu-latest, macos-latest${NC}"
else
    echo -e "${YELLOW}⚠ OS matrix may be different${NC}"
fi

# Extract compilers from CI workflow
COMPILERS=$(grep -A 5 "matrix:" "${WORKFLOWS_DIR}/ci.yml" | grep "compiler:" | sed 's/.*\[\(.*\)\].*/\1/')
echo -e "${GREEN}✓ Compilers: ${COMPILERS}${NC}"

# Extract MPI implementations from CI workflow
MPI=$(grep -A 5 "matrix:" "${WORKFLOWS_DIR}/ci.yml" | grep "mpi:" | sed 's/.*\[\(.*\)\].*/\1/')
echo -e "${GREEN}✓ MPI implementations: ${MPI}${NC}"

echo ""
echo "=== Checking Project Structure ==="

# Check for CMakeLists.txt
if [ -f "${PROJECT_ROOT}/CMakeLists.txt" ]; then
    echo -e "${GREEN}✓ CMakeLists.txt found${NC}"
else
    echo -e "${YELLOW}⚠ CMakeLists.txt not found (using Makefile)${NC}"
fi

# Check for source files
CPP_COUNT=$(find "${PROJECT_ROOT}" -name "*.cpp" -type f | wc -l | tr -d ' ')
H_COUNT=$(find "${PROJECT_ROOT}" -name "*.h" -o -name "*.hpp" | wc -l | tr -d ' ')
echo -e "${GREEN}✓ Found ${CPP_COUNT} .cpp files${NC}"
echo -e "${GREEN}✓ Found ${H_COUNT} .h/.hpp files${NC}"

# Check for tests
if [ -d "${PROJECT_ROOT}/tests" ] || [ -d "${PROJECT_ROOT}/test" ]; then
    echo -e "${GREEN}✓ Tests directory found${NC}"
else
    echo -   -e "${YELLOW}⚠ Tests directory not found${NC}"
fi

echo ""
echo "=== Summary ==="
echo -e "${GREEN}✓ All required workflows are configured${NC}"
echo ""
echo "Next steps:"
echo "1. Commit and push the workflow files"
echo "2. Check the Actions tab in GitHub to see the workflows running"
echo "3. Review the workflow logs for any issues"
echo ""
echo "For more information, see: ${WORKFLOWS_DIR}/README.md"
echo ""
