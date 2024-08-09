from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="foldseek_align",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="A package for running Foldseek analysis and creating distance matrices",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/foldseek_align",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Scientists",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy",
    ],
    extras_require={
        "foldseek": ["foldseek"],
    },
    entry_points={
        "console_scripts": [
            "foldseek_align=foldseek_align.foldseek_align:main",
        ],
    },
)