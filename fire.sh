#!/bin/bash

THIS_DIR=$(realpath $(dirname ${0}))

IGA_REPOS_NAME=iga-ads

IGA_ROOT=$(find ~ -type d -name ${IGA_REPOS_NAME})
IGA_EXAMPLES=${IGA_ROOT}/examples
IGA_WILDFIRE=${IGA_EXAMPLES}/wildfire

PROJECT_WILDFIRE_NAME=project_wildfire
PROJECT_WILDFIRE_ROOT=$(find ~ -type d -name ${PROJECT_WILDFIRE_NAME})
PROJECT_WILDFIRE_FILES=${PROJECT_WILDFIRE_ROOT}/wildfire.*

cp ${PROJECT_WILDFIRE_FILES} ${IGA_WILDFIRE}

IGA_BUILD=build

cmake --build ${IGA_ROOT}/${IGA_BUILD}

PROJECT_WILDFIRE_OUT=${PROJECT_WILDFIRE_ROOT}/fire-out
PROJECT_WILDFIRE_GIF=${PROJECT_WILDFIRE_ROOT}/fire-gif

mkdir -p ${PROJECT_WILDFIRE_OUT}
mkdir -p ${PROJECT_WILDFIRE_GIF}

if [ ! -z "$(ls -A ${PROJECT_WILDFIRE_OUT})" ]; then
    rm ${PROJECT_WILDFIRE_OUT}/*.data
fi

if [ ! -z "$(ls -A ${PROJECT_WILDFIRE_GIF})" ]; then
    rm ${PROJECT_WILDFIRE_GIF}/*.gif
fi

# run wildfire executable from {PROJECT_WILDFIRE_OUT} path to store output data there

LOG_TMP=${THIS_DIR}/step_info.tmp
touch ${LOG_TMP}
truncate --size 0 ${LOG_TMP}

cd ${PROJECT_WILDFIRE_OUT}
${IGA_ROOT}/${IGA_BUILD}/examples/wildfire >> ${LOG_TMP}

wait

mapfile -t indices < <(awk '/Step/ {print $2}' ${LOG_TMP})
DATA_INDEX_FIRST=${indices[0]}
DATA_INDEX_SECOND=${indices[1]}
DATA_INDEX_STEP=$((${DATA_INDEX_SECOND} - ${DATA_INDEX_FIRST}))
DATA_INDEX_LAST=${indices[-1]}

rm ${LOG_TMP}

# plot data
GP_HELPER_NAME=gnuplot_helper.gp
GP_HELPER_PATH=$(find ~ -type f -name ${GP_HELPER_NAME})

GP_FUEL_GIF=fuel.gif
GP_OUT_GIF=out.gif

cd ${PROJECT_WILDFIRE_OUT}
gnuplot -e "output_gif='${GP_FUEL_GIF}'" -e "input_prefix='fuel'" -e "step=${DATA_INDEX_STEP}" -e "last_index=${DATA_INDEX_LAST}" ${GP_HELPER_PATH}
gnuplot -e "output_gif='${GP_OUT_GIF}'" -e "input_prefix='out'" -e "step=${DATA_INDEX_STEP}" -e "last_index=${DATA_INDEX_LAST}" ${GP_HELPER_PATH}

mv ${GP_FUEL_GIF} ${GP_OUT_GIF} ${PROJECT_WILDFIRE_GIF}
