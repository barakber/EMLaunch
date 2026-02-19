.PHONY: test help clean install notebooks compose-up

# Default target
help:
	@echo "EMLaunch - Electromagnetic Launcher Simulation"
	@echo ""
	@echo "Available targets:"
	@echo "  make notebooks - Start Jupyter notebook server"
	@echo "  make test      - Run the test suite"
	@echo "  make install   - Install package dependencies"
	@echo "  make clean     - Clean up generated files"
	@echo "  make compose-up - Rebuild and start Docker container"
	@echo "  make help      - Show this help message"

# Start Jupyter notebook server (via Julia's IJulia)
notebooks:
	@echo "Precompiling EMLaunch and dependencies..."
	julia --project=. -e 'using Pkg; Pkg.precompile()'
	@echo ""
	@echo "Starting JupyterLab via IJulia..."
	@echo "This ensures the Julia kernel is properly configured"
	@echo "Press Ctrl+C to stop the server"
	@echo ""
	julia --project=. -e 'using IJulia; jupyterlab(dir="notebooks")'

# Run tests
test:
	julia --project=. -e 'using Pkg; Pkg.test()'

# Install dependencies
install:
	julia --project=. -e 'using Pkg; Pkg.instantiate()'

# Rebuild and start Docker container
compose-up:
	docker-compose up --build

# Clean generated files
clean:
	@echo "Cleaning generated files..."
	@find . -type f -name "*.png" -path "./examples/*" -delete
	@find . -type f -name "*.log" -delete
	@find . -type d -name ".ipynb_checkpoints" -exec rm -rf {} + 2>/dev/null || true
	@echo "Clean complete"
