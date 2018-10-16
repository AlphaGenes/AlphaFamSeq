
# List all src files below inside set
SET(AlphaFamSeq_src
        ${SRCAlphaFamSeq}ComputeStat.f90
        ${SRCAlphaFamSeq}specFileModule.f90
	${SRCAlphaFamSeq}AlphaVarCallParallelised.f90
        ${SRCAlphaFamSeq}AlphaFamSeqPlusGeneProb.f90
)

add_definitions(-D COMMIT=${VERSION})
ADD_EXECUTABLE(${AlphaFamSeqEXE} ${AlphaFamSeq_src})


# add external links here
TARGET_LINK_LIBRARIES(${AlphaFamSeqEXE} ${GPLIB})
TARGET_LINK_LIBRARIES(${AlphaFamSeqEXE} ${AHLIB})
INSTALL(TARGETS ${AlphaFamSeqEXE} RUNTIME DESTINATION bin)

