import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="reinvent-scoring-rdkit",
    version="0.0.74",
    author="Girish.C",
    author_email="girishchiluveru09@gmail.com",
    description="Scoring functions for Reinvent",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Girish09-work/reinvent-scoring",
    package_data={"reinvent_scoring": ["scoring/score_components/synthetic_accessibility/fpscores.pkl.gz"]},
    packages=setuptools.find_packages(exclude='unittest_reinvent'),
    classifiers=[
        "Programming Language :: Python :: 3",
       
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    
)
