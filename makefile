CC = gcc
CFLAGS = -O3 -march=native -mtune=native -fopenmp -pipe
LDFLAGS = -lm

TARGET = yard_rects
SRC = yard_rects.c

all: $(TARGET)

$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $(SRC) -o $(TARGET) $(LDFLAGS)

clean:
	rm -f $(TARGET)
