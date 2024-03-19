# Guide for developers

## Changing databases

The easiest way to use alternative databases is to download them manually and
provide them as paths in the configuarion.

If you wish to deploy the module with other databases, you can also fork the
repository and modify the post-deploy scripts for each environement
(`<ENV NAME>.post.deploy.sh`).
Simply modifying the URL vairable should be enough in most cases.

## Data validation

Data validation is performed through Pydantic model validations.
