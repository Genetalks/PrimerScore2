#ifndef TM_CALCULATOR_H
#define TM_CALCULATOR_H

#include <string>
#include <memory>
#include "options.h"
#include "thal.h"
#include "thal_parameters.h"

namespace pt {

class tm_calculator {
public:
    tm_calculator(const pt::options &opt);
    ~tm_calculator();

    void alignment_tm(const std::string &align0, const std::string &seq1, const std::string &seq2, double *tm, double *dG, double *dS, double *dH);

public:
    static const int32_t MAX_LOOP = 30;

private:
    thal_args a;
    thal_results o;
    thal_parameters thermodynamic_parameters;
    double saltCorrection;
};

typedef std::shared_ptr<tm_calculator> tm_calculator_ptr;

}

#endif
