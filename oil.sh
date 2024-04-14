#!/bin/bash

TARGET_PROJECT_DIR=oil
TARGET_PROJECT_NAME=oil2d

THIS_DIR=$(realpath $(dirname ${0}))

IGA_REPOS_NAME=iga-ads

IGA_ROOT=$(find ~ -type d -name ${IGA_REPOS_NAME})
IGA_EXAMPLES=${IGA_ROOT}/examples
IGA_OIL=${IGA_EXAMPLES}/${TARGET_PROJECT_DIR}

PROJECT_OIL_NAME=project_oil2d
PROJECT_OIL_ROOT=$(find ~ -type d -name ${PROJECT_OIL_NAME})
PROJECT_OIL_FILES=${PROJECT_OIL_ROOT}/oil2d.*

cp ${PROJECT_OIL_FILES} ${IGA_OIL}

IGA_BUILD=build

cmake --build ${IGA_ROOT}/${IGA_BUILD}

PROJECT_OIL_OUT=${PROJECT_OIL_ROOT}/oil-out
PROJECT_OIL_GIF=${PROJECT_OIL_ROOT}/oil-gif

mkdir -p ${PROJECT_OIL_OUT}
mkdir -p ${PROJECT_OIL_GIF}

if [ ! -z "$(ls -A ${PROJECT_OIL_OUT})" ]; then
    rm ${PROJECT_OIL_OUT}/*.data
fi

if [ ! -z "$(ls -A ${PROJECT_OIL_GIF})" ]; then
    rm ${PROJECT_OIL_GIF}/*.gif
fi

# run oil2d executable from {PROJECT_WILDFIRE_OUT} path to store output data there

LOG_TMP=${THIS_DIR}/step_info.tmp
touch ${LOG_TMP}
truncate --size 0 ${LOG_TMP}

cd ${PROJECT_OIL_OUT}
${IGA_ROOT}/${IGA_BUILD}/examples/${TARGET_PROJECT_NAME} >> ${LOG_TMP}

wait

mapfile -t indices < <(awk '/Step/ {print $2}' ${LOG_TMP})
DATA_INDEX_FIRST=${indices[0]}
DATA_INDEX_SECOND=${indices[1]}

# trim trailing ','
DATA_INDEX_FIRST=$(echo ${DATA_INDEX_FIRST} | tr -d ',')
DATA_INDEX_SECOND=$(echo ${DATA_INDEX_SECOND} | tr -d ',')

DATA_INDEX_STEP=$((${DATA_INDEX_SECOND} - ${DATA_INDEX_FIRST}))

# plot data
GP_HELPER_NAME=gif_oil.gnu
GP_HELPER_PATH=$(find ~ -type f -name ${GP_HELPER_NAME})

GP_OIL_GIF=oil.gif

cd ${PROJECT_OIL_OUT}
gnuplot -e "output_gif='${GP_OIL_GIF}'" -e "input_prefix='out'" -e "step=100" -e "last_index=10000" ${GP_HELPER_PATH}

mv oil_out.gif ${PROJECT_OIL_GIF}
