ROOTPATH=..
include ${ROOTPATH}/compiler.h
include ${ROOTPATH}/library_defines.h
include ${ROOTPATH}/dependencies.h


conserved_training: ${ROOTPATH}/conserved_training/conserved_training.cpp ${ROOTPATH}/conserved_training/conserved_training.h
	${LINK_IN_EXE_DIR} -I/home/yfu7/Applications/ /home/yfu7/Applications/dlib/all/source.cpp -lpthread -lX11 ${ROOTPATH}/conserved_training/conserved_training.cpp ${CONFIG_FILE_PARSER} ${CMD_LINE_PARSER} ${RNA_FILES} ${PHMM_FILES} 

conserved_training_smp: ${ROOTPATH}/conserved_training/conserved_training.cpp ${ROOTPATH}/conserved_training/conserved_training.h
	${LINK_IN_EXE_DIR} -fopenmp -I/home/yfu7/Applications/ /home/yfu7/Applications/dlib/all/source.cpp -lpthread -lX11 ${ROOTPATH}/conserved_training/conserved_training.cpp ${CONFIG_FILE_PARSER} ${CMD_LINE_PARSER} ${RNA_FILES_SMP} ${PHMM_FILES}

