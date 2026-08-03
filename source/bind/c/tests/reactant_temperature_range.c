#include "cea.h"

#include <math.h>

int main(void)
{
    cea_real temperature_min;
    cea_real temperature_max;

    if (cea_reactant_get_valid_temperature_range(
            "(CH2)x(cr)", &temperature_min, &temperature_max) != CEA_INVALID_SIZE) {
        return 1;
    }
    if (cea_init() != CEA_SUCCESS) {
        return 2;
    }
    if (cea_reactant_get_valid_temperature_range(
            "(CH2)x(cr)", &temperature_min, &temperature_max) != CEA_SUCCESS) {
        return 3;
    }
    if (fabs(temperature_min - 288.15) > 1e-12 || fabs(temperature_max - 308.15) > 1e-12) {
        return 4;
    }
    if (cea_reactant_get_valid_temperature_range(
            "not-a-reactant", &temperature_min, &temperature_max) != CEA_INVALID_INDEX) {
        return 5;
    }

    return 0;
}
