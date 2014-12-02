CC=g++
CFLAG=-std=c++11 -O3 -g -fopenmp
metis=/work/01802/mooncake/opt/metis
parmetis=/work/01802/mooncake/opt/parmetis
liblinear=/work/01802/mooncake/opt/liblinear-1.94
#objects = $(addprefix $(liblinear)/, linear.o tron.o)
#objects = $(addprefix $(liblinear)/blas/, *.o)

col:
	$(CC) $(CFLAG) -I$(liblinear) -c collaborative_ranking.cpp 
	$(CC) $(CFLAG) -o a.out collaborative_ranking.o $(liblinear)/linear.o $(liblinear)/tron.o $(liblinear)/blas/*.o
alt:
	$(CC) $(CFLAG) -o alt_rankSVM alt_rankSVM.cpp 
ocs:
	$(CC) $(CFLAG) -I$(metis)/include -I$(parmetis)/include -L$(metis)/lib -L$(parmetis)/lib -o omega2csr omega2csr.cpp -lparmetis -lmetis
run:
	./rate2partition ../dataset/ml-10m.dat movielens m1.g 4
clean:
	rm *.o[0-9]* *.e[0-9]* 
