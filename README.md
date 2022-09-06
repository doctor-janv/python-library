# Doctor-JanV's python-library
My Python library

## How to use:

All python-scripts needing this library should add the following to the top of
the script:
```python
import sys
import os 

def CheckModules(module_name : str):
    module_path = os.getenv(module_name)
    if module_path == None: sys.exit(module_name + " env var undefined!")
    else: sys.path.append(module_path)

CheckModules("DOCTOR_JANV_PYTHON_LIBRARY")
```

The environment variable `DOCTOR_JANV_PYTHON_LIBRARY` should point to the path 
of the library (i.e., the one that **contains** `ScientificComputing`, not within 
`ScientificComputing`).