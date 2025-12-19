# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Run Commands
- Run simulations: `python exp.py`
- Run specific mode: Edit `modes` list in exp.py and run `python exp.py`

## Code Style Guidelines
- Imports: Standard library first, then third-party packages, then local modules
- Type hints: Use Python typing module and include annotations for function parameters and returns
- Variable naming: Use snake_case for variables and functions
- Error handling: Use try/except blocks for error-prone operations
- Use NumPy functions over loops whenever possible for performance
- Document functions with proper docstrings using the NumPy style format
- Use f-strings for string formatting
- Keep line length under 100 characters
- Maintain 4-space indentation
- For numerical simulations, always use fixed random seeds for reproducibility