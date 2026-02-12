![Panl Logo](https://github.com/user-attachments/assets/394c4f2e-60f8-4c6b-8fb5-2df8f39c5de5)

# Panl

[![Tests](https://github.com/Czarified/panl/actions/workflows/tests.yml/badge.svg)](https://github.com/Czarified/panl/actions/workflows/tests.yml)
[![Python Version](https://img.shields.io/badge/python-3.11%2B-blue)](https://www.python.org/downloads/)
[![Poetry](https://img.shields.io/badge/dependency%20management-poetry-blue)](https://python-poetry.org/)

**Panl** is a Python library for analyzing structural panels with cutouts, providing tools for stress analysis, geometric optimization, and design validation.

## Features

- 🔧 Structural panel analysis with cutout support
- 📊 Stress concentration calculations
- 🎯 Geometric optimization
- ✅ Comprehensive testing with pytest
- 🎨 Code quality enforced with pre-commit hooks

## Installation

### From PyPI (once published)

```bash
pip install panl
```

### For Development

```bash
# Clone the repository
git clone https://github.com/Czarified/panl.git
cd panl

# Install Poetry if you haven't already
# Visit: https://python-poetry.org/docs/#installation

# Install dependencies
poetry install

# Install pre-commit hooks
poetry run pre-commit install
```

## Quick Start

```python
import panl

# Your code here
print(f"Panl version: {panl.__version__}")
```

## Development

This project uses several tools to maintain code quality:

- **Poetry**: Dependency management
- **Nox**: Task automation and testing across Python versions
- **pre-commit**: Git hooks for code quality
- **pytest**: Testing framework with coverage

### Running Tests

```bash
# Run tests with coverage
poetry run nox -s tests

# Run tests for specific Python version
poetry run nox -s tests-3.11
```

### Code Quality

```bash
# Run all linters
poetry run nox -s lint

# Format code with black
poetry run nox -s black

# Sort imports with isort
poetry run nox -s isort

# Run pre-commit hooks
poetry run nox -s pre-commit
```

### Available Nox Sessions

```bash
# List all available sessions
poetry run nox --list

# Run a specific session
poetry run nox -s <session-name>
```

## Requirements

- Python >= 3.11
- Poetry for dependency management

## Contributing

Contributions are welcome! Please ensure:

1. All tests pass: `poetry run nox -s tests`
2. Code is formatted: `poetry run nox -s black`
3. Linting passes: `poetry run nox -s lint`
4. Pre-commit hooks pass: `poetry run pre-commit run --all-files`

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Author

Benjamin Crews - your.email@example.com

This project is an Agentic Experiment by the author, and uses Gemini 3 Flash.
