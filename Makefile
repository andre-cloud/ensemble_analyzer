docs: 
	cd docs/
	make html

format:
	/opt/homebrew/Caskroom/miniconda/base/envs/main/bin/python -m autopep8 -r --in-place src
	/opt/homebrew/Caskroom/miniconda/base/envs/main/bin/python -m black src

check:
	/opt/homebrew/Caskroom/miniconda/base/envs/main/bin/python -m autopep8 -rdv src

test:
	/opt/homebrew/Caskroom/miniconda/base/envs/main/bin/python -m coverage run --source=. -m pytest
	/opt/homebrew/Caskroom/miniconda/base/envs/main/bin/python -m coverage html
	# /opt/homebrew/Caskroom/miniconda/base/envs/main/bin/python -m coverage report
	rm *xyz