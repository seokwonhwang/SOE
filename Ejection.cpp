#include "Ejection.h"

using namespace SOENS;

std::ostream& operator<<(std::ostream& os, const Ejection& ejection)
{
    os << "!------------------------ Ejection --------------------------!" << std::endl;
    os << "M_i = " << ejection.mi_ << std::endl;
    os << "M_t = " << ejection.mt_ << std::endl;
    os << "!-------------------- End of Ejection -----------------------!" << std::endl;

    return os;
}
