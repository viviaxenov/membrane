UTest: utask.o
	g++ test/UTest.cpp obj/mem.o obj/utask.o -std=c++11 -O2 -o UTest 2> errors.err
test : mem.o vec.o
	g++ test/test.cpp obj/mem.o obj/vec.o -std=c++11 -O2 -o test/test 2> errors.err
utask.o : mem.o
	[ ! -d "./obj" ] && mkdir obj || :
	g++ source/uniform.cpp obj/mem.o  -c -O2 -o obj/utask.o 2> errors.err
mem.o : 
	[ ! -d "./obj" ] && mkdir obj || :
	g++ source/membrane.cpp -c -O2  -o obj/mem.o 2> errors.err
vec.o : 
	[ ! -d "./obj" ] && mkdir obj || :
	g++ source/vec3.cpp -c -O2  -o obj/vec.o 2> errors.err
clean : 
	[ -d "./obj" ] && rm -r obj || :
