




include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules



convect: convect.o  chkopts
		-${FLINKER} -o convect convect.o  ${PETSC_LIB}
		${RM} convect.o

clean::
	${RM} convect



