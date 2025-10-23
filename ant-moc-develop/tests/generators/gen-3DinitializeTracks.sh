#!/bin/sh

testdir=../autotests
mkdir ${testdir} &> /dev/null

## Templates used to generate tests
base_dir=../unittests
base_src=test_3D_initializeTracks.cpp
base_input=test_3D_initializeTracks_input.inc

base_stem="${base_src%.*}"
testsuf="${base_src#*_}"
teststem="autotest${testsuf%.*}"

for num_azim in 4 8; do
for azim_spacing in 0.5 0.1; do
for num_polar in 4 8; do
for polar_spacing in 0.75 0.5; do
for xmax in 20 32.13; do
for ymax in 20 32.13; do
testsuite="${teststem}INPUTa${num_azim}as${azim_spacing//.}p${num_polar}ps${polar_spacing//.}"
testsuite="${testsuite}x${xmax//./_}y${ymax//./_}"
testfile="${testsuite}.cpp"
testinput="${testsuite}_input.inc"
## Generate tests
cat ${base_dir}/${base_src} > ${testdir}/${testfile}
sed -i "s/${base_stem}/${testsuite}/g" ${testdir}/${testfile}
sed -i "s/${base_input}/${testinput}/g" ${testdir}/${testfile}
## Generate test input
cat << END  > ${testdir}/${testinput}
/* Auto-generated Input */
num_azim = ${num_azim};
azim_spacing = ${azim_spacing};
num_polar = ${num_polar};
polar_spacing = ${polar_spacing};
tolerance = 1e-5;
max_iters = 1;
axial_refines = 5;

xmin = -${xmax};
xmax =  ${xmax};
ymin = -${ymax};
ymax =  ${ymax};
zmin = -32.13;
zmax = 32.13;
/* Auto-generated Input */
END
done
done
done
done
done
done
