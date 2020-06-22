from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='optimalTAD',
      version='0.1',
      url='https://github.com/cosmoskaluga/optimalTAD',
      author='Khrameeva Lab',
      author_email='dmitrii.smirnov@phystech.edu',
      long_description=long_description,
      long_description_content_type="text/markdown",
      license='MIT',
      packages=find_packages(),
      zip_safe=False,
      platforms='any',
      include_package_data=True,
      install_requires=['numpy>=1.16.2',
                        'scipy>=1.2.1',
                        'h5py>=2.9.0'],
      classifiers=[
                   'Development Status :: 3 - Alpha',
                   'License :: OSI Approved :: MIT License',
                   'Programming Language :: Python :: 3.8'
                   ]
      )
