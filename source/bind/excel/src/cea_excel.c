#include "cea_excel.h"

#include <stdio.h>
#include <string.h>

int cea_excel_test_add(int a, int b)
{
    return a + b;
}

int cea_excel_version(char *version, int version_len)
{
    char version_text[64];
    int version_text_len;

    if (version == NULL) {
        return CEA_EXCEL_ERROR_NULL_BUFFER;
    }
    if (version_len <= 0) {
        return CEA_EXCEL_ERROR_INVALID_LENGTH;
    }

    version[0] = '\0';

    version_text_len = snprintf(
        version_text,
        sizeof(version_text),
        "%d.%d.%d",
        CEA_EXCEL_VERSION_MAJOR,
        CEA_EXCEL_VERSION_MINOR,
        CEA_EXCEL_VERSION_PATCH
    );
    if (version_text_len < 0) {
        return CEA_EXCEL_ERROR_TRUNCATED;
    }

    if (version_text_len >= version_len) {
        memcpy(version, version_text, (size_t)(version_len - 1));
        version[version_len - 1] = '\0';
        return CEA_EXCEL_ERROR_TRUNCATED;
    }

    memcpy(version, version_text, (size_t)version_text_len + 1U);
    return CEA_EXCEL_SUCCESS;
}
