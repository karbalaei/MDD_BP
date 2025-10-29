#!/bin/bash

echo "**** Job starts ****"
date

find ./Results/Amygdala/ -type f -name "*_isoTWAS.RDS" | sed 's#./Results/Amygdala/##g; s#_isoTWAS.RDS##g' > Amygdala_RDS_list
find ./Results/Amygdala/ -type f -name "*_isoTWAS.RDS" | sed 's#./Results/Amygdala/##g; s#_isoTWAS.RDS##g' > sACC_RDS_list

sbatch 6_running_burden_test_Amygdala.sh
sbatch 6_running_burden_test_Amygdala_BP.sh
sbatch 6_running_burden_test_sACC.sh
sbatch 6_running_burden_test_sACC_BP.sh

echo "**** Job ends ****"
date