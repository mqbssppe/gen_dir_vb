CXX = g++

CXXFLAGS = -fopenmp

PROGRAMS = gd_new_format gd_log_format dir_new_format dir_log_format gen_dir_rpkm chib


all: $(PROGRAMS)

get_k.o: get_k.h get_k.cpp infiles.h

transposeFiles.o: transposeFiles.cpp transposeFiles.h

common.o: common.h common.cpp

gd_new_format: get_k.o play2.cpp infiles.h
	$(CXX) $(CXXFLAGS) -o gd_new_format get_k.o play2.cpp

gd_log_format: get_k.o logformat.cpp infiles.h
	$(CXX) $(CXXFLAGS) -o gd_log_format get_k.o logformat.cpp

dir_new_format: get_k.o play2_dir.cpp infiles.h
	$(CXX) $(CXXFLAGS) -o dir_new_format get_k.o play2_dir.cpp

dir_log_format: get_k.o logformat_dir.cpp infiles.h
	$(CXX) $(CXXFLAGS) -o dir_log_format get_k.o logformat_dir.cpp


chib: get_k.o gibbs.cpp infiles.h 
	$(CXX) $(CXXFLAGS) -o chib get_k.o gibbs.cpp

gen_dir_rpkm: common.o transposeFiles.o gen_dir_rpkm.cpp get_k.cpp
	$(CXX) $(CXXFLAGS) common.o transposeFiles.o -o gen_dir_rpkm gen_dir_rpkm.cpp get_k.cpp

clean:
	rm *.o $(PROGRAMS)


