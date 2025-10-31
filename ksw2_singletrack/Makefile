CC = gcc
CXX = g++
CFLAGS = -g -Wall -Wextra -O2 -march=native -std=c++17
CPPFLAGS = -DHAVE_KALLOC
INCLUDES = -I. -Iksw2
KSW2_DIR = ksw2
KSW2_OBJS = $(KSW2_DIR)/ksw2_gg.o $(KSW2_DIR)/ksw2_gg2.o $(KSW2_DIR)/ksw2_gg2_sse.o \
            $(KSW2_DIR)/ksw2_extz.o $(KSW2_DIR)/ksw2_extz2_sse.o \
            $(KSW2_DIR)/ksw2_extd.o $(KSW2_DIR)/ksw2_extd2_sse.o \
            $(KSW2_DIR)/ksw2_extf2_sse.o $(KSW2_DIR)/ksw2_exts2_sse.o \
			$(KSW2_DIR)/kalloc.o

OBJS = benchmark_singletrack.o ksw2_singletrack.o
LIBS = -lz
TARGET = benchmark_singletrack

all: $(TARGET)

$(TARGET): $(OBJS) $(KSW2_OBJS)
	$(CXX) $(CFLAGS) $(INCLUDES) $(CPPFLAGS) -o $@ $^ $(LIBS)

benchmark_singletrack.o: benchmark_singletrack.cpp ksw2_singletrack.hpp
	$(CXX) $(CFLAGS) $(INCLUDES) $(CPPFLAGS) -c $< -o $@

ksw2_singletrack.o: ksw2_singletrack.cpp ksw2_singletrack.hpp
	$(CXX) $(CFLAGS) $(INCLUDES) $(CPPFLAGS) -c $< -o $@

$(KSW2_OBJS):
	$(MAKE) -C $(KSW2_DIR)

clean:
	rm -f $(OBJS) $(TARGET)
	$(MAKE) -C $(KSW2_DIR) clean
