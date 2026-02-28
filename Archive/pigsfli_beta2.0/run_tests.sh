#!/bin/bash

# ----------------------------------------
#  Pretty colors
# ----------------------------------------
GREEN="\033[1;32m"
RED="\033[1;31m"
YELLOW="\033[1;33m"
BLUE="\033[1;34m"
RESET="\033[0m"

echo -e "${BLUE}Compiling and running all tests...${RESET}"

# ----------------------------------------
#  Track results
# ----------------------------------------
PASSED=0
FAILED=0
TOTAL=0

# ----------------------------------------
#  Loop over all test_*.cpp files
# ----------------------------------------
for test in tests/test_*.cpp; do
    if [[ -f "$test" ]]; then
        TOTAL=$((TOTAL + 1))
        exe="${test%.cpp}"
        exe="${exe##*/}"

        echo -e "${YELLOW}----------------------------------------${RESET}"
        echo -e "${BLUE}Running $exe${RESET}"

        # Compile
        if g++ -std=c++17 "$test" -I./src -o "$exe"; then
            # Run
            if ./"$exe"; then
                echo -e "${GREEN}[PASS] $exe${RESET}"
                PASSED=$((PASSED + 1))
            else
                echo -e "${RED}[FAIL] $exe (runtime error)${RESET}"
                FAILED=$((FAILED + 1))
            fi
        else
            echo -e "${RED}[FAIL] $exe (compile error)${RESET}"
            FAILED=$((FAILED + 1))
        fi
    fi
done

echo -e "${YELLOW}----------------------------------------${RESET}"
echo -e "${BLUE}Test Summary:${RESET}"
echo -e "  ${GREEN}Passed: $PASSED${RESET}"
echo -e "  ${RED}Failed: $FAILED${RESET}"
echo -e "  ${YELLOW}Total:  $TOTAL${RESET}"

if [[ $FAILED -eq 0 ]]; then
    echo -e "${GREEN}All tests passed successfully.${RESET}"
else
    echo -e "${RED}Some tests failed. Check output above.${RESET}"
fi
