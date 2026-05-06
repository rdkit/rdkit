#include <iostream>
#include <cstdio>

int main() {
    printf("Program started\n");
    fflush(stdout);

    #include "TemplateSmiles.h"

    printf("TemplateSmiles.h included successfully\n");
    printf("Number of templates: %zu\n", TEMPLATE_SMILES.size());
    fflush(stdout);

    return 0;
}
