docs: 
	cd docs/
	make html

format:
	/Users/andrea/opt/miniconda3/bin/python -m autopep8 -r --in-place src
	/Users/andrea/opt/miniconda3/bin/python -m black src

check:
	/Users/andrea/opt/miniconda3/bin/python -m autopep8 -rdv src

test:
	/Users/andrea/opt/miniconda3/bin/python -m coverage run --source=. -m pytest
	/Users/andrea/opt/miniconda3/bin/python -m coverage html
	# /Users/andrea/opt/miniconda3/bin/python -m coverage report
	rm *xyz