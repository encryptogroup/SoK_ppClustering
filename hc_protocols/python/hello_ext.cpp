#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/python.hpp>
#include <string>

struct World
{
    World(std::string msg): msg(msg) {} // added constructor
    void set(std::string msg) { this->msg = msg; }
    std::string greet() { return msg; }
    std::string msg;
};

using namespace boost::python;

BOOST_PYTHON_MODULE(hello_ext)
{
    class_<World>("World", init<std::string>())
        .def("greet", &World::greet)
        .def("set", &World::set)
    ;
}
