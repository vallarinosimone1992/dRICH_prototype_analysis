#!/bin/bash

cmake --build .

OUTDIR="/home/simone/Documents/PhD/Attivita/DRichPrototype/Results/2022_test_beam_passports/aerogel_studies_analysis"

sleep 2

#for I in 66 67 68 69 70 71 72 74 75 76 81 82 83 84 85 120 121 122 123 124 125 130 131 132 133 147 148 149 150 151 152 153 154 155 156 157 158 159 160 161
#for I in 74 75 76 81 82 83 84 85 154 155 
for I in 274 276 277 279 280 281 313 314 315 316
do
  ./reco $I
  echo #
done

:'
NAME="aerogel_JAP4_sideA"
./mon $NAME 66
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP4_sideB"
./mon $NAME 67
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP3_sideA"
./mon $NAME 68
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP3_sideB"
./mon $NAME 69
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP2_sideA"
./mon $NAME 70
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done 
echo #

NAME="aerogel_JAP2_sideB"
./mon $NAME 71
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP2_sideA_JAP3_sideB"
./mon $NAME 72
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP2_sideA_JAP3_sideB_JAP3_sideA"
./mon $NAME 74 75 76 86
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #


NAME="aerogel_JAP2_sideA_JAP3_sideB_JAP3_sideA_20GeV"
./mon $NAME 81 82 85
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP2_sideA_JAP3_sideB_JAP3_sideA_50GeV"
./mon $NAME 83 84
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="no_aerogel"
./mon $NAME 120
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP2_sideA_JAP3_sideB_JAP3_sideA_additionalLucite"
./mon $NAME 121
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP2_sideA_JAP3_sideB_JAP3_sideA_filter_window_330_370"
./mon $NAME 122
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP2_sideA_JAP3_sideB_JAP3_sideA_filter_shortpass_400_sideA"
./mon $NAME 123
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP2_sideA_JAP3_sideB_JAP3_sideA_filter_shortpass_400_sideB"
./mon $NAME 124
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP2_sideA_JAP3_sideB_JAP3_sideA_filter_longpass_400_sideB"
./mon $NAME 125
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="standard_conditions"
./mon $NAME 130 131 132 133
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="standard_conditions_wavelength_filter_window_400_450"
./mon $NAME 147 148 149
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="standard_conditions_wavelength_filter_window_450_500"
./mon $NAME 150 151
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="standard_conditions_wavelength_filter_window_500_550"
./mon $NAME 152 153
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="standard_conditions_wavelength_filter_window_550_600"
./mon $NAME 155
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="standard_conditions_wavelength_filter_longpass_600"
./mon $NAME 156 157
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="standard_conditions_wavelength_filter_window_280_320"
./mon $NAME 158 159
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_JAP1_sideA"
./mon $NAME 160 161
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #
'

NAME="aerogel_RUS10_sideA_10GeV"
./mon $NAME 274 315
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_RUS10_sideB_10GeV"
./mon $NAME 316
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_RUS10_sideA_lucite_300_10GeV"
./mon $NAME 276 277 279 313
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_RUS10_sideA_lucite_300_10GeV"
./mon $NAME 314
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_RUS12_sideA_lucite_300_8GeV"
./mon $NAME 280
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

NAME="aerogel_RUS12_sideA_lucite_300_4GeV"
./mon $NAME 281
mv /home/simone/Work/EIC/dRICH/prototype/testBeam//output/plot/${NAME}/passport.pdf ${OUTDIR}/passport_${NAME}.pdf
echo $NAME done
echo #

