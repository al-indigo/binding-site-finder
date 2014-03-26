#include "functions.h"

void init_ahoc(AhoCorasickPlus& atm, Pwm& matrix, size_t patterns_allowed, size_t& total_words, optimization_type type) {
    std::vector<std::vector<char> > patterns;
    
    patterns = matrix.getWords(patterns_allowed, type);
    
    for (size_t i = 0; i < patterns.size(); i++)
    {
        AhoCorasickPlus::EnumReturnStatus res;
        AhoCorasickPlus::PatternId patId = total_words++;
        res = atm.addPattern(patterns[i], patId);
        if (res!=AhoCorasickPlus::RETURNSTATUS_SUCCESS)
        {
            std::cerr << "Failed to add: " << std::endl;
        }
    }

    patterns.clear();
}