# dRICH prototype reco and analysis
Repository which includes all the corrections and analysis for the dRICH prototype.

* * *

## How to set the environment
1. Set the environment variable $DRICH_SUITE
2. Clone this repository inside $DRICH_SUITE
3. Run the setenv.sh script
4. In $DRICH_SUITE/DATA create a software link (ln -s) between dRICH_DATA and the path to the directory which contain the prototype data. (Or simply, put directly the dRICH_DATA directory inside the $DRICH_SUITE/DATA);
5. In $DRICH_SUITE/DATA create a software link (ln -s) between GEM_DATA and the path to the directory which contain the tracking data. (Again as point 4)
6. In $DRICH_SUITE/DATA create a directory named logbook that will include the file logbook.tsv.
7. Then you are ready to compile the software.

* * *

## How to compile
1. cd $DRICH_SUITE/dRICH_prototype_analysis/sw/build
2. rm -r *
3. cmake ..
4. cmake --build .
5. This produce two executable, reco and ana.

* * *

## How to run
1. Be sure that $DRICH_SUITE/DATA/logbook/logbook.tsv exists and contains info of specific run you want to analyze. Maybe you have to download it and rename it froma google sheets; in this case, remember to download it in tsv format.
2. cd $DRICH_SUITE/dRICH_prototype_analysis/sw/build
3. ./reco followed by the run number. In any case, I hope I had time to write the help option.
4. The reconstruction produces the ROOT TTree with integrated data in $DRICH_SUITE/processed_data/integrated_dRICH_GEM_data/
5. ./ana followed by: the output file name, one or more run number. You can also run with just one argument, the run number. In this case, you will produce the plots only for the selected run, in a directory with the run number name.
