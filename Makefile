install:
	RUSTFLAGS="-Awarnings" cargo build  --release
	mkdir bin
	mv target/release/voxbirch ./bin
	mv target/release/vbirchlow ./bin

clean:
	cargo clean && rm -rf ./bin

fresh:
	make clean && make install
