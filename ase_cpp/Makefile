EXES=src/Ase

DEPS=

release: deps/installed.txt
	make -j24 -C build/$@
	cd build/$@ && cp $(EXES) ../../bin

debug: deps/installed.txt
	make -j24 -C build/$@

clean:
	make -C build/release clean
	make -C build/debug clean
	rm -rf bin/*

$(DEPS):
	mkdir -p deps/$@/build
	cd deps/$@/build && cmake .. -DCMAKE_INSTALL_PREFIX=../../../src && make && make install && cd ../../.. 
	# Now anything in deps/$@/include or lib will now be in src

deps/installed.txt: $(DEPS)
	touch deps/installed.txt
