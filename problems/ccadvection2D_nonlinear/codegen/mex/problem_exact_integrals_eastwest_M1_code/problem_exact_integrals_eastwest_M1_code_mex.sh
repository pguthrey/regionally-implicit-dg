MATLAB="/Applications/MATLAB_R2018a.app"
Arch=maci64
ENTRYPOINT=mexFunction
MAPFILE=$ENTRYPOINT'.map'
PREFDIR="/Users/piersonguthrey/Library/Application Support/MathWorks/MATLAB/R2018a"
OPTSFILE_NAME="./setEnv.sh"
. $OPTSFILE_NAME
COMPILER=$CC
. $OPTSFILE_NAME
echo "# Make settings for problem_exact_integrals_eastwest_M1_code" > problem_exact_integrals_eastwest_M1_code_mex.mki
echo "CC=$CC" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "CFLAGS=$CFLAGS" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "CLIBS=$CLIBS" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "COPTIMFLAGS=$COPTIMFLAGS" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "CDEBUGFLAGS=$CDEBUGFLAGS" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "CXX=$CXX" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "CXXFLAGS=$CXXFLAGS" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "CXXLIBS=$CXXLIBS" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "CXXOPTIMFLAGS=$CXXOPTIMFLAGS" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "CXXDEBUGFLAGS=$CXXDEBUGFLAGS" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "LDFLAGS=$LDFLAGS" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "LDOPTIMFLAGS=$LDOPTIMFLAGS" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "LDDEBUGFLAGS=$LDDEBUGFLAGS" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "Arch=$Arch" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "LD=$LD" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo OMPFLAGS= >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo OMPLINKFLAGS= >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "EMC_COMPILER=clang" >> problem_exact_integrals_eastwest_M1_code_mex.mki
echo "EMC_CONFIG=optim" >> problem_exact_integrals_eastwest_M1_code_mex.mki
"/Applications/MATLAB_R2018a.app/bin/maci64/gmake" -j 1 -B -f problem_exact_integrals_eastwest_M1_code_mex.mk
