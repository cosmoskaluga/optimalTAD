from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='optimalTAD',
      version='0.1.0',
      url='https://github.com/cosmoskaluga/optimalTAD',
      author='Khrameeva Lab',
      author_email='Dmitrii.Smirnov@skoltech.ru',
      long_description=long_description,
      long_description_content_type="text/markdown",
      entry_points={'console_scripts': ['optimalTAD = optimalTAD.__main__:optimalTAD']},
      license='MIT',
      packages=find_packages(),
      include_package_data=True,
      zip_safe=False,
      platforms='any',
      install_requires=['numpy>=1.16.2',
                        'scipy>=1.2.1',
                        'pandas>=1.0.3',
                        'h5py>=2.9.0',
                        'seaborn>=0.10.0',
                        'pyBigWig'],
      classifiers=[
                   'Development Status :: 3 - Alpha',
                   'License :: OSI Approved :: MIT License',
                   'Programming Language :: Python :: 3.8', 
                   'Programming Language :: Python :: 3.9']
      )
