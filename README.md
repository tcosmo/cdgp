# Welcome!

# Commands
## Git related commands

- When work has been done on a file: `git add file` or `git add -A` (will add all modified/untracked files)
- When work needs to be recorded: `git commit -m "short message"`
- To push to the repo: `git push`
- To get work from others: `git pull`

## Python related commands

### Virtual env
- In order to activate the virtual environment: `source venv/bin/activate`
- To quit the environment: either quit terminal, or type `deactivate`

### Documentation
- `cd docs`
- `make html` produces the html doc in folder `build/`
- The result is in `docs/build/html/index.html`

If you want to run doc tests:
- `make doctest` runs all documentation tests



### Tests
- `python -m unittest -q tests.example_tests`