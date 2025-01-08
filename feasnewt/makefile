
KNITRO_HOME = /home/sbrantq

METRIC = imr

.SUFFIXES: .c .o

#CC = icc
#CFLAGS = -O3 -DREORDER -DCHECK_MESH
#LFLAGS =

#CC = pgcc
#CFLAGS = -O4
#LFLAGS = 

# Original configuration
# CC = gcc
# #CFLAGS = -O9 -Wall -DCHECK_MESH -DSEGMENT_MESH -DREORDER -DSWEEP
# CFLAGS = -O9 -Wall -fPIC -DCHECK_MESH -DREORDER -DCHECK_PD  -DSWEEP
# LFLAGS = -lm

CC = clang
CFLAGS = -O3 -Wall -fPIC -DCHECK_MESH -DREORDER -DCHECK_PD  -DSWEEP -ffast-math -fno-finite-math-only
LFLAGS = -lm

#CC = cc
#CFLAGS = -xO5 -dalign 
#LFLAGS = -lm

HSRC  = main/meshOpt.c
HSRC += elements/${METRIC}/fcn3i_am.c elements/elem3.c
HSRC += meshes/hex.c meshes/mesh3.c meshes/mesh.c
HSRC += optimizers/pre3.c 
HSRC += generator/gen3.c 

HLSRC += elements/${METRIC}/fcn3i_am.c elements/elem3.c
HLSRC += meshes/hex.c meshes/mesh3.c meshes/mesh.c
HLSRC += optimizers/pre3.c 
HLSRC += generator/gen3.c 

LHSRC  = main/meshOpt.c
LHSRC += elements/${METRIC}/fcn3i_am.c elements/${METRIC}/fcn3il_am.c 
LHSRC += elements/elem3.c elements/elem3l.c
LHSRC += meshes/hex.c meshes/mesh3.c meshes/mesh.c
LHSRC += optimizers/opt3_gs.c
LHSRC += generator/gen3.c 

QHSRC  = main/meshOpt.c
QHSRC += elements/${METRIC}/fcn3id_am.c elements/elem3d.c
QHSRC += meshes/hex.c meshes/mesh3.c meshes/mesh.c
QHSRC += optimizers/opt3_qn.c
QHSRC += optimizers/pre3.c 
QHSRC += generator/gen3.c 

LHLSRC += elements/${METRIC}/fcn3i_am.c elements/${METRIC}/fcn3il_am.c 
LHLSRC += elements/elem3.c elements/elem3l.c
LHLSRC += meshes/hex.c meshes/mesh3.c meshes/mesh.c
LHLSRC += optimizers/opt3_gs.c
LHLSRC += generator/gen3.c 

HOSRC  = main/meshConvert.c
HOSRC += elements/${METRIC}/fcn3i_am.c elements/${METRIC}/fcn3i_gm.c 
HOSRC += elements/${METRIC}/fcn3i_hm.c elements/elem3.c
HOSRC += meshes/hex.c meshes/mesh3.c meshes/mesh.c

TSRC  = main/meshOpt.c
TSRC += elements/${METRIC}/fcn3e_am.c elements/elem3.c
TSRC += meshes/tet.c meshes/mesh3.c meshes/mesh.c
TSRC += optimizers/pre3.c 
TSRC += generator/gen3.c 

TLSRC += elements/${METRIC}/fcn3e_am.c elements/elem3.c
TLSRC += meshes/tet.c meshes/mesh3.c meshes/mesh.c
TLSRC += optimizers/pre3.c 
TLSRC += generator/gen3.c 

LTSRC  = main/meshOpt.c
LTSRC += elements/${METRIC}/fcn3e_am.c elements/${METRIC}/fcn3el_am.c 
LTSRC += elements/elem3.c elements/elem3l.c
LTSRC += meshes/tet.c meshes/mesh3.c meshes/mesh.c
LTSRC += optimizers/opt3_gs.c
LTSRC += generator/gen3.c 

QTSRC  = main/meshOpt.c
QTSRC += elements/${METRIC}/fcn3ed_am.c elements/elem3d.c
QTSRC += meshes/tet.c meshes/mesh3.c meshes/mesh.c
QTSRC += optimizers/opt3_qn.c
QTSRC += optimizers/pre3.c 
QTSRC += generator/gen3.c 

LTLSRC += elements/${METRIC}/fcn3e_am.c elements/${METRIC}/fcn3el_am.c 
LTLSRC += elements/elem3.c elements/elem3l.c
LTLSRC += meshes/tet.c meshes/mesh3.c meshes/mesh.c
LTLSRC += optimizers/opt3_gs.c
LTLSRC += generator/gen3.c 

TCSRC  = main/meshCheck.c
TCSRC += elements/${METRIC}/fcn3e_am.c elements/${METRIC}/fcn3e_gm.c 
TCSRC += elements/${METRIC}/fcn3e_hm.c 
TCSRC += elements/cond/fcn3e_am.c elements/elem3.c
TCSRC += meshes/tet.c meshes/mesh3.c meshes/mesh.c
TCSRC += check/check.c

TOSRC  = main/meshConvert.c
TOSRC += elements/${METRIC}/fcn3e_am.c elements/${METRIC}/fcn3e_gm.c 
TOSRC += elements/${METRIC}/fcn3e_hm.c elements/elem3.c
TOSRC += meshes/tet.c meshes/mesh3.c meshes/mesh.c

QSRC  = main/meshOpt.c
QSRC += elements/${METRIC}/fcn2i_am.c elements/elem2.c
QSRC += meshes/quad.c meshes/mesh2.c meshes/mesh.c
QSRC += optimizers/pre2.c 

LQSRC  = main/meshOpt.c
LQSRC += elements/${METRIC}/fcn2i_am.c elements/${METRIC}/fcn2il_am.c 
LQSRC += elements/elem2.c elements/elem2l.c
LQSRC += meshes/quad.c meshes/mesh2.c meshes/mesh.c
LQSRC += optimizers/opt2_gs.c

QOSRC  = main/meshConvert.c
QOSRC += elements/${METRIC}/fcn2i_am.c elements/${METRIC}/fcn2i_gm.c 
QOSRC += elements/${METRIC}/fcn2i_hm.c elements/elem2.c
QOSRC += meshes/quad.c meshes/mesh2.c meshes/mesh.c

RSRC  = main/meshOpt.c
RSRC += elements/${METRIC}/fcn2e_am.c elements/elem2.c
RSRC += meshes/tri.c meshes/mesh2.c meshes/mesh.c
RSRC += optimizers/pre2.c 

LRSRC  = main/meshOpt.c
LRSRC += elements/${METRIC}/fcn2e_am.c elements/${METRIC}/fcn2el_am.c 
LRSRC += elements/elem2.c elements/elem2l.c
LRSRC += meshes/tri.c meshes/mesh2.c meshes/mesh.c
LRSRC += optimizers/opt2_gs.c

RCSRC  = main/meshCheck.c
RCSRC += elements/${METRIC}/fcn2e_am.c elements/${METRIC}/fcn2e_gm.c 
RCSRC += elements/${METRIC}/fcn2e_hm.c 
RCSRC += elements/cond/fcn2e_am.c elements/elem2.c
RCSRC += meshes/tri.c meshes/mesh2.c meshes/mesh.c
RCSRC += check/check.c

ROSRC  = main/meshConvert.c
ROSRC += elements/${METRIC}/fcn2e_am.c elements/${METRIC}/fcn2el_am.c 
ROSRC += elements/elem2.c elements/elem2l.c
ROSRC += meshes/tri.c meshes/mesh2.c meshes/mesh.c

# Turn on the CHECK_PD flag when the metric chose is not guaranteed to
# be convex.

ifeq (${METRIC},dft)
	CFLAGS += -DCHECK_PD
endif

ifeq (${METRIC},dftt)
	CFLAGS += -DCHECK_PD
endif
 
all: opt check refine convert lib

opt: tetOpt hexOpt triOpt quadOpt
check: tetCheck triCheck
refine: tetRefine hexRefine triRefine quadRefine
convert: tetConvert hexConvert triConvert quadConvert
knitro: knitro2 knitro3
lib: tetLib hexLib

tetOpt: tetOpt_ls tetOpt_tr tetOpt_trls tetOpt_gs tetOpt_qn
hexOpt: hexOpt_ls hexOpt_tr hexOpt_trls hexOpt_gs hexOpt_qn
triOpt: triOpt_ls triOpt_tr triOpt_trls triOpt_gs
quadOpt: quadOpt_ls quadOpt_tr quadOpt_trls quadOpt_gs
knitro2: triOpt_ktr triOpt_ktrh quadOpt_ktr quadOpt_ktrh
knitro3: tetOpt_ktr tetOpt_ktrh hexOpt_ktr hexOpt_ktrh

tetLib: libtetOpt_ls.so libtetOpt_tr.so libtetOpt_trls.so libtetOpt_gs.so
hexLib: libhexOpt_ls.so libhexOpt_tr.so libhexOpt_trls.so libhexOpt_gs.so

tetOpt_ls: optimizers/opt3_ls.c ${TSRC}
	${CC} ${CFLAGS} -Iinclude -o tetOpt_ls optimizers/opt3_ls.c ${TSRC} ${LFLAGS}

tetOpt_tr: optimizers/opt3_tr.c ${TSRC}
	${CC} ${CFLAGS} -Iinclude -o tetOpt_tr optimizers/opt3_tr.c ${TSRC} ${LFLAGS}

tetOpt_trls: optimizers/opt3_trls.c ${TSRC}
	${CC} ${CFLAGS} -Iinclude -o tetOpt_trls optimizers/opt3_trls.c ${TSRC} ${LFLAGS}

tetOpt_gs: ${LTSRC}
	${CC} ${CFLAGS} -Iinclude -o tetOpt_gs ${LTSRC} ${LFLAGS}

tetOpt_qn: ${QTSRC}
	${CC} ${CFLAGS} -Iinclude -o tetOpt_qn ${QTSRC} ${LFLAGS}

libtetOpt_ls.so: optimizers/opt3_ls.c ${TLSRC}
#	cat ${TLSRC} optimizers/opt3_ls.c > tetOpt_ls.c
#	${CC} ${CFLAGS} -Iinclude -o tetOpt_ls tetOpt_ls.c ${LFLAGS}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -shared -o libtetOpt_ls.so optimizers/opt3_ls.c ${TLSRC} ${LFLAGS}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -c optimizers/opt3_ls.c ${TLSRC} 
	ar ruv libtetOpt_ls.a *.o
	/bin/rm -f *.o

libtetOpt_tr.so: optimizers/opt3_tr.c ${TLSRC}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -shared -o libtetOpt_tr.so optimizers/opt3_tr.c ${TLSRC} ${LFLAGS}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -c optimizers/opt3_tr.c ${TLSRC} 
	ar ruv libtetOpt_tr.a *.o
	/bin/rm -f *.o

libtetOpt_trls.so: optimizers/opt3_trls.c ${TLSRC}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -shared -o libtetOpt_trls.so optimizers/opt3_trls.c ${TLSRC} ${LFLAGS}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -c optimizers/opt3_trls.c ${TLSRC} 
	ar ruv libtetOpt_trls.a *.o
	/bin/rm -f *.o

libtetOpt_gs.so: ${LTLSRC}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -shared -o libtetOpt_gs.so ${LTLSRC} ${LFLAGS}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -c ${LTLSRC} 
	ar ruv libtetOpt_gs.a *.o
	/bin/rm -f *.o

tetOpt_ktr: optimizers/opt3_ktr.c ${TSRC}
	${CC} ${CFLAGS} -Iinclude -I${KNITRO_HOME}/include -o tetOpt_ktr optimizers/opt3_ktr.c ${TSRC} -L${KNITRO_HOME}/lib -lknitro ${LFLAGS} -ldl

tetOpt_ktrh: optimizers/opt3_ktrh.c ${TSRC}
	${CC} ${CFLAGS} -Iinclude -I${KNITRO_HOME}/include -o tetOpt_ktrh optimizers/opt3_ktrh.c ${TSRC} -L${KNITRO_HOME}/lib -lknitro ${LFLAGS} -ldl

tetCheck: ${TCSRC}
	${CC} ${CFLAGS} -DCHECK -Iinclude -o tetCheck ${TCSRC} ${LFLAGS}

tetRefine: ${TSRC}
	${CC} ${CFLAGS} -DREFINE -Iinclude -o tetRefine ${TSRC} ${LFLAGS}

tetConvert: ${TOSRC}
	${CC} ${CFLAGS} -Iinclude -o tetConvert ${TOSRC} ${LFLAGS}

hexOpt_ls: optimizers/opt3_ls.c ${HSRC}
	${CC} ${CFLAGS} -Iinclude -o hexOpt_ls optimizers/opt3_ls.c ${HSRC} ${LFLAGS}

hexOpt_tr: optimizers/opt3_tr.c ${HSRC}
	${CC} ${CFLAGS} -Iinclude -o hexOpt_tr optimizers/opt3_tr.c ${HSRC} ${LFLAGS}

hexOpt_trls: optimizers/opt3_trls.c ${HSRC}
	${CC} ${CFLAGS} -Iinclude -o hexOpt_trls optimizers/opt3_trls.c ${HSRC} ${LFLAGS}

hexOpt_gs: ${LHSRC}
	${CC} ${CFLAGS} -Iinclude -o hexOpt_gs ${LHSRC} ${LFLAGS}

hexOpt_qn: ${QHSRC}
	${CC} ${CFLAGS} -Iinclude -o hexOpt_qn ${QHSRC} ${LFLAGS}

libhexOpt_ls.so: optimizers/opt3_ls.c ${HLSRC}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -shared -o libhexOpt_ls.so optimizers/opt3_ls.c ${HLSRC} ${LFLAGS}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -c optimizers/opt3_ls.c ${HLSRC} 
	ar ruv libhexOpt_ls.a *.o
	/bin/rm -f *.o

libhexOpt_tr.so: optimizers/opt3_tr.c ${HLSRC}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -shared -o libhexOpt_tr.so optimizers/opt3_tr.c ${HLSRC} ${LFLAGS}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -c optimizers/opt3_tr.c ${HLSRC} 
	ar ruv libhexOpt_tr.a *.o
	/bin/rm -f *.o

libhexOpt_trls.so: optimizers/opt3_trls.c ${HLSRC}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -shared -o libhexOpt_trls.so optimizers/opt3_trls.c ${HLSRC} ${LFLAGS}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -c optimizers/opt3_trls.c ${HLSRC} 
	ar ruv libhexOpt_trls.a *.o
	/bin/rm -f *.o

libhexOpt_gs.so: ${LHLSRC}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -shared -o libhexOpt_gs.so ${LHLSRC} ${LFLAGS}
	${CC} ${CFLAGS} -DLIBRARY -Iinclude -c ${LHLSRC} 
	ar ruv libhexOpt_gs.a *.o
	/bin/rm -f *.o

hexOpt_ktr: optimizers/opt3_ktr.c ${HSRC}
	${CC} ${CFLAGS} -Iinclude -I${KNITRO_HOME}/include -o hexOpt_ktr optimizers/opt3_ktr.c ${HSRC} -L${KNITRO_HOME}/lib -lknitro ${LFLAGS} -ldl

hexOpt_ktrh: optimizers/opt3_ktrh.c ${HSRC}
	${CC} ${CFLAGS} -Iinclude -I${KNITRO_HOME}/include -o hexOpt_ktrh optimizers/opt3_ktrh.c ${HSRC} -L${KNITRO_HOME}/lib -lknitro ${LFLAGS} -ldl

hexRefine: ${HSRC}
	${CC} ${CFLAGS} -DREFINE -Iinclude -o hexRefine ${HSRC} ${LFLAGS}

hexConvert: ${HOSRC}
	${CC} ${CFLAGS} -Iinclude -o hexConvert ${HOSRC} ${LFLAGS}

triOpt_ls: optimizers/opt2_ls.c ${RSRC}
	${CC} ${CFLAGS} -Iinclude -o triOpt_ls optimizers/opt2_ls.c ${RSRC} ${LFLAGS}

triOpt_tr: optimizers/opt2_tr.c ${RSRC}
	${CC} ${CFLAGS} -Iinclude -o triOpt_tr optimizers/opt2_tr.c ${RSRC} ${LFLAGS}

triOpt_trls: optimizers/opt2_trls.c ${RSRC}
	${CC} ${CFLAGS} -Iinclude -o triOpt_trls optimizers/opt2_trls.c ${RSRC} ${LFLAGS}

triOpt_gs: ${LRSRC}
	${CC} ${CFLAGS} -Iinclude -o triOpt_gs ${LRSRC} ${LFLAGS}

triOpt_ktr: optimizers/opt2_ktr.c ${RSRC}
	${CC} ${CFLAGS} -Iinclude -I${KNITRO_HOME}/include -o triOpt_ktr optimizers/opt2_ktr.c ${RSRC} -L${KNITRO_HOME}/lib -lknitro ${LFLAGS} -ldl

triOpt_ktrh: optimizers/opt2_ktrh.c ${RSRC}
	${CC} ${CFLAGS} -Iinclude -I${KNITRO_HOME}/include -o triOpt_ktrh optimizers/opt2_ktrh.c ${RSRC} -L${KNITRO_HOME}/lib -lknitro ${LFLAGS} -ldl

triCheck: ${RCSRC}
	${CC} ${CFLAGS} -DCHECK -Iinclude -o triCheck ${RCSRC} ${LFLAGS}

triRefine: ${RSRC}
	${CC} ${CFLAGS} -DREFINE -Iinclude -o triRefine ${RSRC} ${LFLAGS}

triConvert: ${ROSRC}
	${CC} ${CFLAGS} -Iinclude -o triConvert ${ROSRC} ${LFLAGS}

quadOpt_ls: optimizers/opt2_ls.c ${QSRC}
	${CC} ${CFLAGS} -Iinclude -o quadOpt_ls optimizers/opt2_ls.c ${QSRC} ${LFLAGS}

quadOpt_tr: optimizers/opt2_tr.c ${QSRC}
	${CC} ${CFLAGS} -Iinclude -o quadOpt_tr optimizers/opt2_tr.c ${QSRC} ${LFLAGS}

quadOpt_trls: optimizers/opt2_trls.c ${QSRC}
	${CC} ${CFLAGS} -Iinclude -o quadOpt_trls optimizers/opt2_trls.c ${QSRC} ${LFLAGS}

quadOpt_gs: ${LQSRC}
	${CC} ${CFLAGS} -Iinclude -o quadOpt_gs ${LQSRC} ${LFLAGS}

quadOpt_ktr: optimizers/opt2_ktr.c ${QSRC}
	${CC} ${CFLAGS} -Iinclude -I${KNITRO_HOME}/include -o quadOpt_ktr optimizers/opt2_ktr.c ${QSRC} -L${KNITRO_HOME}/lib -lknitro ${LFLAGS} -ldl

quadOpt_ktrh: optimizers/opt2_ktrh.c ${QSRC}
	${CC} ${CFLAGS} -Iinclude -I${KNITRO_HOME}/include -o quadOpt_ktrh optimizers/opt2_ktrh.c ${QSRC} -L${KNITRO_HOME}/lib -lknitro ${LFLAGS} -ldl

quadConvert: ${QOSRC}
	${CC} ${CFLAGS} -Iinclude -o quadConvert ${QOSRC} ${LFLAGS}

quadRefine: ${QSRC}
	${CC} ${CFLAGS} -DREFINE -Iinclude -o quadRefine ${QSRC} ${LFLAGS}

clean:
	/bin/rm -f tetOpt_ktr tetOpt_ktrh tetOpt_ls tetOpt_tr tetOpt_trls tetOpt_gs tetOpt_qn tetCheck tetConvert
	/bin/rm -f hexOpt_ktr hexOpt_ktrh hexOpt_ls hexOpt_tr hexOpt_trls hexOpt_gs hexOpt_qn hexConvert
	/bin/rm -f triOpt_ktr triOpt_ktrh triOpt_ls triOpt_tr triOpt_trls triOpt_gs triCheck triConvert
	/bin/rm -f quadOpt_ktr quadOpt_ktrh quadOpt_ls quadOpt_tr quadOpt_trls quadOpt_gs quadConvert
	/bin/rm -f triRefine tetRefine quadRefine hexRefine
	/bin/rm -f lib*.so lib*.a
	/bin/rm -f core core.* *.o *.il *~

