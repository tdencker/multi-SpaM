all: copy

.PHONY: copy clean
copy:
	make -C ./src
	cp ./src/multi-SpaM ./bin
clean:
	make -C ./src clean