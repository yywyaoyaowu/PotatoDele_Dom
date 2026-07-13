output=$1
vcf=$2
bed=$3

SIMULATE="ReLERNN_SIMULATE"
TRAIN="ReLERNN_TRAIN"
PREDICT="ReLERNN_PREDICT"
BSCORRECT="ReLERNN_BSCORRECT"
SEED="42"
MU="1e-8"
URTR="10"
DIR="./${output}"
VCF="./${vcf}"
GENOME="./${bed}"

source /home/huyong/software/anaconda3/bin/activate ML_python3.9
module load cuda/cuda-12.2

# Simulate data
${SIMULATE} \
    --vcf ${VCF} \
    --genome ${GENOME} \
    --projectDir ${DIR} \
    --assumedMu ${MU} \
    --upperRhoThetaRatio ${URTR} \
    --nTrain 13000 \
    --forceWinSize 15000 \
	--nVali 2000 \
    --nTest 100 \
    --nCPU 3 \
    --seed ${SEED}

# Train network
${TRAIN} \
    --projectDir ${DIR} \
    --nCPU 3 \
    --seed ${SEED}

# Predict
${PREDICT} \
    --vcf ${VCF} \
    --projectDir ${DIR} \
    --seed ${SEED}
