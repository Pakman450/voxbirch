install:
	RUSTFLAGS="-Awarnings" cargo build  --release
	mkdir bin
	mv target/release/voxbirch ./bin

clean:
	cargo clean && rm -rf ./bin

fresh:
	make clean && make install

test:
	RUST_TEST_THREADS=1 cargo test

deny:
	cargo deny check licenses

about:
	cargo about generate about.hbs --config about.toml -o THIRD_PARTY_NOTICES.html
