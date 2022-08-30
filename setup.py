from setuptools import find_packages, setup

setup(
    name="q2-humann3",
    version="0.0.0",
    packages=find_packages(),
    # install_requires=[
    #     "qiime >= 2.0.0",
    #     "humann3 >= 0.9.4, < 1.0.0",
    #     "biom-format >= 2.1.5, < 2.2.0",
    # ],
    author="",
    author_email="",
    description="QIIME2 plugin for running HUMAnN3",
    entry_points={"qiime2.plugins": ["q2-humann3=q2_humann3.plugin_setup:plugin"]},
)
