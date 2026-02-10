
help:
	@echo "make clean   - Clean the build artifacts and remove the ./bin directory"
	@echo "make deny    - Check for license issues using cargo-deny"
	@echo "make about   - Generate THIRD_PARTY_NOTICES.html using cargo-about"

clean:
	cargo clean && rm -rf ./bin

deny:
	cargo deny check licenses

about:
	cargo about generate about.hbs --config about.toml -o THIRD_PARTY_NOTICES.html
