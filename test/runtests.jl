using SecondQuantizedAlgebra
using ParallelTestRunner: ParallelTestRunner

# Start with autodiscovered tests
testsuite = ParallelTestRunner.find_tests(@__DIR__)

# Parse arguments
args = ParallelTestRunner.parse_args(ARGS)

if ParallelTestRunner.filter_tests!(testsuite, args)
    delete!(testsuite, "quality/JET")
    delete!(testsuite, "quality/JET_unitary")
end

ParallelTestRunner.runtests(SecondQuantizedAlgebra, args; testsuite)
