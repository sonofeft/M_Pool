

python setup.py sdist

rem only upload source to PyPI
rem C:\Users\Charlie\AppData\Roaming\Python\Python37\Scripts\twine.exe upload --skip-existing dist/*.gz
C:\Users\Charlie\AppData\Roaming\Python\Scripts\twine.exe upload --skip-existing dist/*.gz
