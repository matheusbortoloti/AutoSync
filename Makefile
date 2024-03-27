.PHONY: build test clean auto_sync parse generate 

#build: setup
#	gcc main.c AutoSync.c -o Main.o -lpthread

parse: 
	@echo "C File to be parsed:" $(C_FILE)	
	@echo "Starting parser...\n"
	python3 src/parser_auto_sync.py $(C_FILE)
	@echo "Parser generated JSON file\n"

generate_code: 
	@echo "Starting AutoSync code generator...\n"
	python3 src/code_generator_auto_sync.py $(C_FILE)

run_c_code:
	@echo "Running generated C code...\n"
	cd 05_Workspace && make all	
	cd ..

auto_sync: parse generate

test: 
	./Main.o	

clean:
	rm -f 05_Workspace/*.c
	rm -f 05_Workspace/*.h
	rm -f 05_Workspace/*.o

	unlink pycparser/examples/parser.py
	unlink pycparser/examples/code_generator