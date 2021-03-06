from setuptools import setup

setup(name='pythonQE',
      version='2.0',
      description='A simple python interface for Quantum Espresso',
      url='https://github.com/ftherrien/pythonQE',
      author='Felix Therrien',
      author_email='felix.therrien@gmail.com',
      license='MIT',
      packages=['pythonQE'],
      install_requires=[
        'numpy',
        'pylada'
        ],
      zip_safe=False)
