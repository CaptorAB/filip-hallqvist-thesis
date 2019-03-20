#ifndef __EMSCRIPTEN__
#define CATCH_CONFIG_RUNNER
#include <lib/catch.h>
#endif

#include <lib/random.h>
#include <src/genetic.cpp>

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
using namespace emscripten;

EMSCRIPTEN_BINDINGS(libcapgen)
{
  value_object<OptimizeOptions>("OptimizeOptions")
      .field("population", &OptimizeOptions::population)
      .field("elitism", &OptimizeOptions::elitism)
      .field("generations", &OptimizeOptions::generations)
      .field("bits", &OptimizeOptions::bits)
      .field("steps", &OptimizeOptions::steps)
      .field("instruments", &OptimizeOptions::instruments);

  value_object<Result>("Result")
      .field("fitness", &Result::fitness)
      .field("chromosome", &Result::chromosome)
      .field("individual", &Result::individual);

  emscripten::register_vector<int>("VectorInt");
  emscripten::register_vector<double>("VectorDouble");

  function("optimize", &optimize);
}

#else

using Random = effolkronium::random_static;
int main(int argc, char *argv[])
{
  Random::seed(42);
  int result = Catch::Session().run(argc, argv);
  return result;
}

#endif