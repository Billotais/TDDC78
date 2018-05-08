

sequential:main_sequential.c  physics.c
	gcc -std=c99 -g -o $@ $^ -lrt -lm
parallel: main_parallel.cpp physics.c
	mpicxx  -std=c++11 -trace -g -o $@ $^ -lrt -lm -Nmpi -tcollect -I$(VT_ROOT)/include -L$(VT_LIB_DIR) -lVT $(VT_ADD_LIBS)

parallel_home: main_parallel.cpp physics.c
	mpicxx -std=c++11 -g -o $@ $^ -lrt -lm
allclean:
	-rm *.o sequential
