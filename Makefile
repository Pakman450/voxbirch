
help:
	@echo "make install - Build the project and move the binary to ./bin"
	@echo "make bindings - Build the project and install the Python bindings"
	@echo "make clean   - Clean the build artifacts and remove the ./bin directory"
	@echo "make fresh   - Clean and then install and build binary"
	@echo "make test    - Run tests with a single thread"
	@echo "make deny    - Check for license issues using cargo-deny"
	@echo "make about   - Generate THIRD_PARTY_NOTICES.html using cargo-about"

install:
	RUSTFLAGS="-Awarnings" cargo build  --release
	mkdir bin
	mv target/release/voxbirch ./bin

bindings:
	cd python && \
	python3 -m venv .venv && \
	source .venv/bin/activate && \
	pip install maturin && \
	maturin develop

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
