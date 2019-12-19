import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

# fetch version from the main code file
__version__ = ''
with open('rim/__init__.py') as initpy:
  for line in initpy:
      if line.startswith('__version__'):
          exec(line.strip())
          break

setuptools.setup(
     name='rim',
     version=__version__,
     author="Ivan Majic",
     author_email="imajicos@gmail.com",
     description="A package for calculating RIM between three spatial objects",
     long_description=long_description,
     long_description_content_type="text/markdown",
     url="https://github.com/SelfHealingMapsProject/RIM",
     packages=setuptools.find_packages(),
     install_requires=['numpy>=1.17.4', 'Shapely>=1.6.4'],
     classifiers=[
         "Programming Language :: Python :: 3",
         "License :: OSI Approved :: MIT License",
         "Operating System :: OS Independent",
     ],
    python_requires='>=3.6',
 )
