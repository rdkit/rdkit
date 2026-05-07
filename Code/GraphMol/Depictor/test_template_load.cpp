#include <iostream>
#include "Templates.h"

int main() {
    std::cout << "Loading templates..." << std::endl;
    auto& templates = RDDepict::CoordinateTemplates::getRingSystemTemplates();
    std::cout << "Templates loaded successfully!" << std::endl;

    // Try to access a template
    if (templates.hasTemplateOfSize(6)) {
        std::cout << "Has templates of size 6" << std::endl;
        auto& t = templates.getMatchingTemplates(6);
        std::cout << "Got " << t.size() << " templates of size 6" << std::endl;
    }

    std::cout << "Test completed!" << std::endl;
    return 0;
}
