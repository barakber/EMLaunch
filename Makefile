.PHONY: test help clean install notebooks

# Default target
help:
	@echo "EMLaunch - Electromagnetic Launcher Simulation"
	@echo ""
	@echo "Available targets:"
	@echo "  make notebooks - Start Jupyter notebook server"
	@echo "  make test      - Run the test suite"
	@echo "  make install   - Install package dependencies"
	@echo "  make clean     - Clean up generated files"
	@echo "  make help      - Show this help message"

# Start Jupyter notebook server (via Julia's IJulia)
notebooks:
	@echo "Starting Jupyter notebook server via IJulia..."
	@echo "This ensures the Julia kernel is properly configured"
	@echo "Press Ctrl+C to stop the server"
	@echo ""
	julia --project=. -e 'using IJulia; notebook(dir="notebooks")'

# Run tests
test:
	julia --project=. -e 'using Pkg; Pkg.test()'

# Install dependencies
install:
	julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Clean generated files
clean:
	@echo "Cleaning generated files..."
	@find . -type f -name "*.png" -path "./examples/*" -delete
	@find . -type f -name "*.log" -delete
	@find . -type d -name ".ipynb_checkpoints" -exec rm -rf {} + 2>/dev/null || true
	@echo "Clean complete"
