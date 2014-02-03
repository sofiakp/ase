DEPS=

default: release

release debug: deps/installed.txt
	make -j4 -C build/$@

clean:
	make -j4 -C build/release clean
	make -j4 -C build/debug clean

$(DEPS):
	mkdir -p deps/$@/build
	cd deps/$@/build && cmake .. -DCMAKE_INSTALL_PREFIX=../../../src && make && make install && cd ../../.. 
	# Now anything in deps/$@/include or lib will now be in src

deps/installed.txt: $(DEPS)
	touch deps/installed.txt

install: release
	make install -C build/release
